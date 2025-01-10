////////////////////////////////////////////////////////////////////////////////
// Copyright (C) 2018 Olivier Delaneau, University of Lausanne
// 
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
// 
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
// 
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.
////////////////////////////////////////////////////////////////////////////////
#include <containers/haplotype_set.h>

haplotype_set::haplotype_set() {
	clear();
}

haplotype_set::~haplotype_set() {
	clear();
}

void haplotype_set::clear() {
	n_site = 0;
	n_hap = 0;
}

void haplotype_set::readGrouping(string ifile) {
	vrb.title("Reading target groups");
	string buffer;
	vector < string > tokens, tokensG, tokensS;
	map < string, int > :: iterator it;
	input_file fd(ifile);
	if (fd.fail()) vrb.error("Cannot open [" + ifile + "]");
	//getline(fd, buffer, '\n');
	while (getline(fd, buffer, '\n')) {
		if (stb.split(buffer, tokens) >= 2) {
			it = IDs_map.find(tokens[0]);
			if (it == IDs_map.end()) vrb.bullet("Sample [" + tokens[0] + "] unfound !");
			else {
				targetIDXs.push_back(it->second);
				sourceIDXs.push_back(vector < vector < int > > ());
				groupIDs.push_back(vector < string > ());

				for (int g = 1 ; g < tokens.size() ; g++) {
					stb.split(tokens[g], tokensG, "=");
					if (tokensG.size() != 2) vrb.error("Could not extract a group ID for [" + tokens[0] + "]");
					groupIDs.back().push_back(tokensG[0]);
					stb.split(tokensG[1], tokensS, ";");
					sourceIDXs.back().push_back(vector < int > ());
					for (int t = 0 ; t < tokensS.size() ; t ++) {
						it = IDs_map.find(tokensS[t]);
						if (it == IDs_map.end()) vrb.bullet("Sample [" + tokens[0] + "] has an undefined source sample [" + tokensS[t] + "] being unfound !");
						else sourceIDXs.back().back().push_back(it->second);
					}
				}
			}
		}
	}
	fd.close();

	hasHole = (IDs_map.size()*2 != n_hap);
	copyingProbabilities = vector < vector < float > > (2*targetIDXs.size());
	for (int t =0 ; t < targetIDXs.size() ; t ++) {
		copyingProbabilities[2*t+0] = vector < float > ((hasHole + sourceIDXs[t].size())*n_site*2, 0.0);
		copyingProbabilities[2*t+1] = vector < float > ((hasHole + sourceIDXs[t].size())*n_site*2, 0.0);
	}
}

void haplotype_set::writeCopyingProbabilities(string ofile, variant_map & V) {
	output_file fd (ofile);
	vrb.bullet("Writing Copying Probabilities in [" + ofile + "]");

	//Write file header
	fd << "#CHROM\tPOS\tIDX\tCM";
	for (int t = 0 ; t < sourceIDXs.size() ; t ++) {
		//Hap0
		for (int g = 0 ; g < sourceIDXs[t].size() ; g ++) fd << "\t" << IDs[targetIDXs[t]] << "_" << groupIDs[t][g] << "_0";
		if (hasHole) fd << "\t" << IDs[targetIDXs[t]] << "_HOLE_0";
		//Hap1
		for (int g = 0 ; g < sourceIDXs[t].size() ; g ++) fd << "\t" << IDs[targetIDXs[t]] << "_" << groupIDs[t][g] << "_1";
		if (hasHole) fd << "\t" << IDs[targetIDXs[t]] << "_HOLE_1";
		
	}
	fd << endl;

	//Write file body
	for (int l = 0 ; l  < n_site ; l ++) {
		fd << V.vec_pos[l]->chr << "\t" << V.vec_pos[l]->bp << "\t" << l << "\t" << stb.str(V.vec_pos[l]->cm, 3);
		for (int t = 0 ; t < sourceIDXs.size() ; t ++) {
			int ngroups = sourceIDXs[t].size() + hasHole;
			for (int g = 0 ; g < ngroups ; g ++) fd << "\t" << stb.str(copyingProbabilities[2*t+0][ngroups*l+g], 2);
			for (int g = 0 ; g < ngroups ; g ++) fd << "\t" << stb.str(copyingProbabilities[2*t+1][ngroups*l+g], 2);
		}
		fd << endl;
	}
	fd.close();
}

void haplotype_set::writeCopyingProbabilitiesBCF(string ofile, variant_map & V) {
	std::string file_format="w";
	if (ofile.size() > 6 && ofile.substr(ofile.size()-6) == "vcf.gz") file_format = "wz";
	if (ofile.size() > 3 && ofile.substr(ofile.size()-3) == "bcf")  file_format = "wb";

	//--- Open output file ---//
    htsFile * out_fp = hts_open(ofile.c_str(),file_format.c_str());
	if ( out_fp == NULL ) vrb.error("Can't write in [" + ofile + "]");
	vrb.bullet("Writing Copying Probabilities in [" + ofile + "] (BCF format)");

	//--- Prepare and write header ---//
    bcf_hdr_t *out_hdr = bcf_hdr_init("w");
    bcf_hdr_append(out_hdr, "##FORMAT=<ID=DS,Number=1,Type=Float,Description=\"Dosage of alternate allele\">");
	std::string contig_name = "##contig=<ID=" + V.vec_pos[0]->chr + ",length=" + std::to_string(V.vec_pos[V.vec_pos.size()-1]->bp) +">";
    bcf_hdr_append(out_hdr, contig_name.c_str());

	// Add each column
	int col_nbr = 0;
	for (int t = 0 ; t < sourceIDXs.size() ; t ++) {
		//Hap0
		for (int g = 0 ; g < sourceIDXs[t].size() ; g ++) {
			std::string new_sample = IDs[targetIDXs[t]] + "_" + groupIDs[t][g] + "_0";
			bcf_hdr_add_sample(out_hdr, new_sample.c_str());
			col_nbr++;
		}
		if (hasHole){
			std::string new_sample = IDs[targetIDXs[t]] + "_UNR_0";
			bcf_hdr_add_sample(out_hdr, new_sample.c_str());
			col_nbr++;
		}
		for (int g = 0 ; g < sourceIDXs[t].size() ; g ++) {
			std::string new_sample = IDs[targetIDXs[t]] + "_" + groupIDs[t][g] + "_1";
			bcf_hdr_add_sample(out_hdr, new_sample.c_str());
			col_nbr++;
		}
		if (hasHole){
			std::string new_sample = IDs[targetIDXs[t]] + "_UNR_1";
			bcf_hdr_add_sample(out_hdr, new_sample.c_str());
			col_nbr++;
		}
	}
    bcf_hdr_add_sample(out_hdr, NULL);      // Finalize header
    //write header copied before
	if (bcf_hdr_write(out_fp, out_hdr)) vrb.error("Failed to write header to output file");

	//--- Prepare and write records to output file ---//
	bcf1_t *out_rec = bcf_init();
	float ds[col_nbr];

	for (int l = 0 ; l  < n_site ; l ++) {
		//out_rec->rid = stoi(V.vec_pos[l]->chr);
		out_rec->rid = bcf_hdr_name2id(out_hdr, V.vec_pos[l]->chr.c_str());
		out_rec->pos= V.vec_pos[l]->bp -1;
		bcf_update_id(out_hdr, out_rec, V.vec_pos[l]->id.c_str());
		std::string ref_alt = V.vec_pos[l]->ref + "," + V.vec_pos[l]->alt;
		bcf_update_alleles_str(out_hdr, out_rec, ref_alt.c_str()); // REF=T, ALT=C
		int d_pos = 0;
		for (int t = 0 ; t < sourceIDXs.size() ; t ++) {
			int ngroups = sourceIDXs[t].size() + hasHole;
			for (int g = 0 ; g < ngroups ; g ++) ds[d_pos++] = std::stod(stb.str(copyingProbabilities[2*t+0][ngroups*l+g], 2));
			for (int g = 0 ; g < ngroups ; g ++) ds[d_pos++] = std::stod(stb.str(copyingProbabilities[2*t+1][ngroups*l+g], 2));
		}
		bcf_update_format_float(out_hdr, out_rec, "DS", &ds, col_nbr);
		if(bcf_write(out_fp, out_hdr, out_rec)< 0) vrb.error("Failing to write VCF/record");
	}
	//close all the stuffs
	bcf_hdr_destroy(out_hdr);
	bcf_close(out_fp);
	bcf_destroy(out_rec);
}

void haplotype_set::computeIbdProbabilities(variant_map & V){
	vrb.bullet("Finding IBD segments");
	//--- Prepare parentalPhase ---//
	parentalPhase = vector < vector < int >> (sourceIDXs.size());
	for(int i = 0; i<sourceIDXs.size(); i++) parentalPhase[i] = vector<int> (n_site);

	//--- Prepare output ---///
	//Write file header

	//--- Find IBD segments ---//
	// for each individuals
	std::string chr = V.vec_pos[0]->chr;
	char P_letter[] = {'A','B','C','D'};
	for (int t = 0 ; t < sourceIDXs.size() ; t ++) {
		// for each site
		int ibd_status=-1; // 0 = A, 1= B , 2=C, 3=D
		int ngroups = sourceIDXs[t].size() + hasHole;
		double last_cm_pos = V.vec_pos[0]->cm;
		int last_bp_pos = V.vec_pos[0]->bp;
		int bp_pos_before = last_bp_pos;
		float PD_G1_sum = 0.0, PD_G2_sum = 0.0;
		int last_l = 0;
		for (int l = 0 ; l  < n_site ; l ++) {
			double cm_pos = V.vec_pos[l]->cm;
			int bp_pos = V.vec_pos[l]->bp;

			//--- Compute probabilities ---//
			// BE CAREFUL OF CASE WHERE YOU HAVE NA //
			float H0G1 = copyingProbabilities[2*t+0][ngroups*l+0]; 
			float H0G2 = copyingProbabilities[2*t+0][ngroups*l+1]; 
			float H0U = copyingProbabilities[2*t+0][ngroups*l+2]; 
			float H1G1 = copyingProbabilities[2*t+1][ngroups*l+0]; 
			float H1G2 = copyingProbabilities[2*t+1][ngroups*l+1]; 
			float H1U = copyingProbabilities[2*t+1][ngroups*l+2]; 

			float PA =  H0G1 * H1G2 + H0G1 * H1U + H1G2 * H0U; // IBD and phasing good
			float PB =  H0G2 * H1G1 + H0G2 * H1U + H1G1 * H0U; // IBD and phasing error
			float PC = H0U * H1U; // no ibd in region
			float PD = H0G1*H1G1 + H0G2*H1G2; // disomy
			float P_sum = PA + PB + PC + PD;

			// normalize
			PA = PA/P_sum;
			PB = PB/P_sum;
			PC = PC/P_sum;
			PD = PD/P_sum;
			float Probs[] = {PA, PB, PC, PD};

			// get PD sum in case we want UPD
			PD_G1_sum+=H0G1*H1G1;
			PD_G2_sum+=H0G2*H1G2;

			//std::cout << l << " " << bp_pos << " " << cm_pos << " " << H0G1 << " " << H0G2 << " " << H0U << " " << H1G1 << " " << H1G2 << " " << H1U << std::endl;

			//--- Check if we are in a new site (IBD status) --//
			// find index of maximum value 
			int curr_ibd_status = 0;
			float max_val = -1;
			for(int i = 0; i<4; i++){ if(max_val < Probs[i]) {max_val=Probs[i]; curr_ibd_status=i;}};
			//std::cout << max_val << " " << PA << " " << PB << " " << PC << " " << PD << " " << curr_ibd_status << " " << ibd_status << " " << last_bp_pos << " " << bp_pos << " " << ibd_status << " " << cm_pos - last_cm_pos << std::endl;

			// if IBD status change, how long was the segments?
			if (curr_ibd_status!=ibd_status && ibd_status >= 0){
				float length = cm_pos - last_cm_pos;
				// find UPD field
				std::string UPD = "U";
				if (ibd_status == 3) UPD = (PD_G1_sum > PD_G2_sum) ? "G1" : "G2";
				// output
				last_cm_pos = cm_pos; last_bp_pos = bp_pos; 
				// store segment
				// you can put a condition with the segment size if you want to filter small segments
				// previous line by theoule: for(int sub_l = last_l; sub_l<l; sub_l++) parentalPhase[t][sub_l] = ibd_status;

				// new line with cM length filter
				for(int sub_l = last_l; sub_l<l; sub_l++){ if(length>=3) {parentalPhase[t][sub_l] = ibd_status;}};
				// end of new line
			
				PD_G1_sum = 0.0; PD_G2_sum=0.0;
				last_l = l;
			}
			ibd_status = curr_ibd_status;
			bp_pos_before = bp_pos;

		}

		
		// last iteration
		std::string UPD = "U";
		if (ibd_status == 3) UPD = (PD_G1_sum > PD_G2_sum) ? "G1" : "G2";
		for(int sub_l = last_l; sub_l<n_site; sub_l++) parentalPhase[t][sub_l] = ibd_status;

	}


}

//1. compute probabilities A, B and ...
void haplotype_set::computeIbdProbabilities(string ofile, variant_map & V){
	vrb.bullet("Finding IBD segments and writing in [" + ofile + "]");
	//--- Prepare parentalPhase ---//
	parentalPhase = vector < vector < int >> (sourceIDXs.size());
	for(int i = 0; i<sourceIDXs.size(); i++) parentalPhase[i] = vector<int> (n_site);

	//--- Prepare output ---///
	//Write file header
	output_file fd (ofile);
	fd << "CHR\tstart\tend\tProb\tlength_CM\ttarget\tUPD";
	fd << endl;

	//--- Find IBD segments ---//
	// for each individuals
	std::string chr = V.vec_pos[0]->chr;
	char P_letter[] = {'A','B','C','D'};
	for (int t = 0 ; t < sourceIDXs.size() ; t ++) {
		// for each site
		int ibd_status=-1; // 0 = A, 1= B , 2=C, 3=D
		int ngroups = sourceIDXs[t].size() + hasHole;
		double last_cm_pos = V.vec_pos[0]->cm;
		int last_bp_pos = V.vec_pos[0]->bp;
		int bp_pos_before = last_bp_pos;
		float PD_G1_sum = 0.0, PD_G2_sum = 0.0;
		int last_l = 0;
		for (int l = 0 ; l  < n_site ; l ++) {
			double cm_pos = V.vec_pos[l]->cm;
			int bp_pos = V.vec_pos[l]->bp;

			//--- Compute probabilities ---//
			// BE CAREFUL OF CASE WHERE YOU HAVE NA //
			float H0G1 = copyingProbabilities[2*t+0][ngroups*l+0]; 
			float H0G2 = copyingProbabilities[2*t+0][ngroups*l+1]; 
			float H0U = copyingProbabilities[2*t+0][ngroups*l+2]; 
			float H1G1 = copyingProbabilities[2*t+1][ngroups*l+0]; 
			float H1G2 = copyingProbabilities[2*t+1][ngroups*l+1]; 
			float H1U = copyingProbabilities[2*t+1][ngroups*l+2]; 

			float PA =  H0G1 * H1G2 + H0G1 * H1U + H1G2 * H0U; // IBD and phasing good
			float PB =  H0G2 * H1G1 + H0G2 * H1U + H1G1 * H0U; // IBD and phasing error
			float PC = H0U * H1U; // no ibd in region
			float PD = H0G1*H1G1 + H0G2*H1G2; // disomy
			float P_sum = PA + PB + PC + PD;

			// normalize
			PA = PA/P_sum;
			PB = PB/P_sum;
			PC = PC/P_sum;
			PD = PD/P_sum;
			float Probs[] = {PA, PB, PC, PD};

			// get PD sum in case we want UPD
			PD_G1_sum+=H0G1*H1G1;
			PD_G2_sum+=H0G2*H1G2;

			//std::cout << l << " " << bp_pos << " " << cm_pos << " " << H0G1 << " " << H0G2 << " " << H0U << " " << H1G1 << " " << H1G2 << " " << H1U << std::endl;

			//--- Check if we are in a new site (IBD status) --//
			// find index of maximum value 
			int curr_ibd_status = 0;
			float max_val = -1;
			for(int i = 0; i<4; i++){ if(max_val < Probs[i]) {max_val=Probs[i]; curr_ibd_status=i;}};
			//std::cout << max_val << " " << PA << " " << PB << " " << PC << " " << PD << " " << curr_ibd_status << " " << ibd_status << " " << last_bp_pos << " " << bp_pos << " " << ibd_status << " " << cm_pos - last_cm_pos << std::endl;

			// if IBD status change, how long was the segments?
			if (curr_ibd_status!=ibd_status && ibd_status >= 0){
				float length = cm_pos - last_cm_pos;
				// find UPD field
				std::string UPD = "U";
				if (ibd_status == 3) UPD = (PD_G1_sum > PD_G2_sum) ? "G1" : "G2";
				// output
				fd << chr << "\t" << last_bp_pos << "\t" << bp_pos_before << "\t" << P_letter[ibd_status] << "\t" << length << "\t" << IDs[targetIDXs[t]] << "\t" << UPD;
				fd << endl;
				last_cm_pos = cm_pos; last_bp_pos = bp_pos; 
				// store segment
				// you can put a condition with the segment size if you want to filter small segments
				for(int sub_l = last_l; sub_l<l; sub_l++) parentalPhase[t][sub_l] = ibd_status;

				PD_G1_sum = 0.0; PD_G2_sum=0.0;
				last_l = l;
			}
			ibd_status = curr_ibd_status;
			bp_pos_before = bp_pos;

		}

		
		// last iteration
		std::string UPD = "U";
		if (ibd_status == 3) UPD = (PD_G1_sum > PD_G2_sum) ? "G1" : "G2";
		fd << chr << "\t" << last_bp_pos << "\t" << V.vec_pos[V.size()-1]->bp << "\t" << P_letter[ibd_status] << "\t" << V.vec_pos[V.size()-1]->cm - last_cm_pos << "\t" << IDs[targetIDXs[t]] << "\t" << UPD;
		fd << endl;
		for(int sub_l = last_l; sub_l<n_site; sub_l++) parentalPhase[t][sub_l] = ibd_status;

	}
	fd.close();

}

void haplotype_set::writeParentPhasedBcf(string ofile, variant_map & V) {
	std::string file_format="w";
	if (ofile.size() > 6 && ofile.substr(ofile.size()-6) == "vcf.gz") file_format = "wz";
	if (ofile.size() > 3 && ofile.substr(ofile.size()-3) == "bcf")  file_format = "wb";

	//--- Open output file ---//
    htsFile * out_fp = hts_open(ofile.c_str(),file_format.c_str());
	if ( out_fp == NULL ) vrb.error("Can't write in [" + ofile + "]");
	vrb.bullet("Writing IBD-based haplotype scaffold (H0 = G1, H1 = G2) in [" + ofile + "]");

	//--- Prepare and write header ---//
    bcf_hdr_t *out_hdr = bcf_hdr_init("w");
	bcf_hdr_append(out_hdr, "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">");
	std::string contig_name = "##contig=<ID=" + V.vec_pos[0]->chr + ",length=" + std::to_string(V.vec_pos[V.vec_pos.size()-1]->bp) +">";
    bcf_hdr_append(out_hdr, contig_name.c_str());

	// Add each column
	int col_nbr = 0;
	for (int t = 0 ; t < sourceIDXs.size() ; t ++) bcf_hdr_add_sample(out_hdr, IDs[targetIDXs[t]].c_str());      // Finalize header
    bcf_hdr_add_sample(out_hdr, NULL);      // Finalize header
    //write header copied before
	if (bcf_hdr_write(out_fp, out_hdr)) vrb.error("Failed to write header to output file");

	//--- Prepare and write records to output file ---//
	bcf1_t *out_rec = bcf_init();
	int * genotypes_out = (int*)malloc(sourceIDXs.size()*2*sizeof(int));

	for (int l = 0 ; l  < n_site ; l ++) {
		out_rec->rid = bcf_hdr_name2id(out_hdr, V.vec_pos[l]->chr.c_str());
		out_rec->pos= V.vec_pos[l]->bp -1;
		bcf_update_id(out_hdr, out_rec, V.vec_pos[l]->id.c_str());
		std::string ref_alt = V.vec_pos[l]->ref + "," + V.vec_pos[l]->alt;
		bcf_update_alleles_str(out_hdr, out_rec, ref_alt.c_str()); // REF=T, ALT=C
		for (int t = 0 ; t < sourceIDXs.size() ; t ++ ){
			int allele_a = 2 + 2 * (int) H_opt_var.get(l, targetIDXs[t]*2);
			int allele_b = 2 + 2 * (int) H_opt_var.get(l, targetIDXs[t]*2+1);
			if(parentalPhase[t][l] > 1) { genotypes_out[t*2] = allele_a; genotypes_out[t*2+1] = allele_b; } // no phasing
			else if (parentalPhase[t][l] == 0) { genotypes_out[t*2] = bcf_gt_phased(bcf_gt_allele(allele_a)); genotypes_out[t*2+1] = bcf_gt_phased(bcf_gt_allele(allele_b));} // h0 h1
			else if (parentalPhase[t][l] == 1) { genotypes_out[t*2] = bcf_gt_phased(bcf_gt_allele(allele_b)); genotypes_out[t*2+1] = bcf_gt_phased(bcf_gt_allele(allele_a));} // h1 h0
			else{ vrb.error("Parental phasinig unknown"); }
		}
		bcf_update_genotypes(out_hdr, out_rec, genotypes_out, sourceIDXs.size()*2);
		if(bcf_write(out_fp, out_hdr, out_rec)< 0) vrb.error("Failed to write VCF/record");
	}

	//close all the stuffs
	bcf_hdr_destroy(out_hdr);
	bcf_close(out_fp);
	bcf_destroy(out_rec);
}
