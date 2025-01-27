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
#include <thorin/thorin_header.h>

void thorin::write_files_and_finalise() {
	vrb.title("Finalization:");

	//step0: multi-threading
	if (options["thread"].as < int > () > 1) pthread_mutex_destroy(&mutex_workers);

	//step1: writing basic output
	std::string output_probs = options["output"].as <string> ();
	if ((output_probs.size() > 3 && output_probs.substr(output_probs.size()-3) == "bcf") | (output_probs.size() > 6 && output_probs.substr(output_probs.size()-6) == "vcf.gz")) H.writeCopyingProbabilitiesBCF(options["output"].as <string> (), V);
	
	else H.writeCopyingProbabilities(options["output"].as < string > (), V);

	if (options.count("ibd")) H.computeIbdProbabilities(options["ibd"].as <string> (), V, options["scaffold-cM"].as <float>());
		// 2. output vcf
	if (options.count("scaffold")){
		if (!options.count("ibd")) H.computeIbdProbabilities(V);
		H.writeParentPhasedBcf(options["scaffold"].as <string> (), V);
	}

	//step2: Measure overall running time
	vrb.bullet("Total running time = " + stb.str(tac.abs_time()) + " seconds");
}
