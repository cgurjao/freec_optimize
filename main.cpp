/*************************************************************************
Copyright (c) 2010-2011, Valentina BOEVA.

>>> SOURCE LICENSE >>>
This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation (www.fsf.org); either version 2 of the
License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

A copy of the GNU General Public License is available at
http://www.fsf.org/licensing/licenses

>>> END OF LICENSE >>>
*************************************************************************/


// main.cpp

#include "SVfinder.h"
#include "GenomeCopyNumber.h"
#include "version.h"
#include "Minipileup.h"
#include <iomanip>
#include <sstream>

using namespace std ;

int verbose = false;
double minMappabilityPerWindow = 0.85;
bool uniqueMatch = false;

static void print_version()
{
    std::ostringstream ostr;
	ostr << "Control-FREEC v" << std::fixed << std::setprecision(1) << FREEC_VERSION << " : calling copy number alterations and LOH regions using deep-sequencing data\n";
	std::cout << ostr.str();
}

static const char* get_conf_file(int argc, char *argv[])
{
  if (argc < 3) {
	std::cerr << "\n\tPlease specify a config file\n\n";
	usage();
	exit(0);
  }

  if (argc > 3 || (strcmp(argv[1], "-conf") != 0 && strcmp(argv[1], "--conf") != 0)) {
	usage();
	exit(0);
  }

  const char* conf_file = argv[2];
  ifstream ifile(conf_file);

  if (!ifile) {
	std::cerr << "\n\tCould not find your config file.. Please, check the existance of "<< conf_file <<"\n\n";
	exit(-1);
  }

  return conf_file;
}

static void thread_init(unsigned int max_threads, unsigned int thread_verbose)
{
  if (max_threads > 1) {
	std::cout << "MT-mode using " << max_threads << " threads\n";
  } else {
	std::cout << "Non MT-mode\n";
  }

  ThreadPoolManager::init(max_threads, thread_verbose);
}

int main(int argc, char *argv[])
{

    print_version();

	const char* conf_file = get_conf_file(argc, argv);

	ConfigFile cf(conf_file);

	//read parameters and initial variables:

	unsigned int max_threads = (int)cf.Value("general", "maxThreads", 1);
	bool thread_verbose = (bool)cf.Value("general", "threadVerbose", "false");

	thread_init(max_threads, thread_verbose ? 1 : 0);

	int minCNAlength = (int)cf.Value("general","minCNAlength", 1);
	cout << "..Minimal CNA length (in windows) is "<< minCNAlength<< "\n";

    bool makeminipileup = (bool)cf.Value("general", "makeminipileup", "false");
    bool tryOtherPloidy = (bool)cf.Value("general", "tryOtherPloidy", "false");
    cerr << tryOtherPloidy;

    std::string sex = (std::string)cf.Value("general","sex", "");
	if (sex.compare("XX") == 0) {
	  std::cout << "..consider the sample being female\n";
	} else if (sex.compare("XY") == 0) {
	  std::cout << "..consider the sample being male\n";
	} else if (sex.compare("") != 0){
	  std::cerr << "Error: \"sex\" can be either XX or XY\n";
	  return 0;
	}

	double breakPointThreshold = (double)cf.Value("general","breakPointThreshold", .8);
	if (breakPointThreshold < 0) {
	  cerr << "\n\n\t!!ERROR!! (but don't be afraid :)\n\n";
	  cerr << "Starting from FREEC v.4.2 we use the threshold on the slope of the slope of the RSSs (instead of simply slope in FREEC v.<4.1) to define number of breakpoints in segmentation. ";
	  cerr << "This method is more robust and should provide a more uniform segmentation for different chromosomes.\n";
	  cerr << "\n\tWe recomend to use \"breakPointThreshold=0.8\"\n\n";
	  cerr << "It should be a positive value. The higher it is, the less breakpoints you will get.\n";
	  cerr << "\n\tI am sorry, but you need to change this value in your config profile.. Or your can just comment it with #, then the default values of 0.8 will be applied\n";
	  return 0;
	}
	std::cout << "..Breakpoint threshold for segmentation of copy number profiles is "<< breakPointThreshold<< "\n";

	int degree = (int)cf.Value("general", "degree", 3);
	std::cout << "..Polynomial degree for \"ReadCount ~ GC-content\" or \"Sample ReadCount ~ Control ReadCount\" is "<< degree<< "\n";

    int teloCentroFlanks = (int)cf.Value("general","telocentromeric", TELO_CENTRO_FLANCS);
	std::cout << "..telocenromeric set to "<<teloCentroFlanks<<"\n";

    bool ifBedGraphOutPut = (bool)cf.Value("general","BedGraphOutput", "false");
	if (!ifBedGraphOutPut) {
	  std::cout << "..FREEC is not going to output normalized copy number profiles into a BedGraph file (for example, for visualization in the UCSC GB). Use \"[general] BedGraphOutput=TRUE\" if you want a BedGraph file\n";
	}

 	bool contaminationAdjustment = (bool)cf.Value("general","contaminationAdjustment", "false");
	if (!contaminationAdjustment) {
	  std::cout << "..FREEC is not going to adjust profiles for a possible contamination by normal cells\n";
	} else {
      std::cout << "..FREEC is going to adjust profiles for a possible contamination by normal cells\n";
      std::cout <<"..set contaminationAdjustment=FALSE if you don't want to use this option because you think that there is no contamiantion of your tumor sample by normal cells (e.g., it is a cell line, or it non-cancer DNA used without a control sample)\n";
	}
	float knownContamination = (float)cf.Value("general","contamination",0);

    if (knownContamination>0) {
            if (contaminationAdjustment == false) {
                    std::cerr << "..set contaminationAdjustment=TRUE if you want to use \"contamination=...\"\n";
                    exit(0);
            }
            if (knownContamination>100 ) {
                    std::cerr << "..contamination should not be greater than 100%\n";
                    exit(0);
            }
            if (knownContamination>1 ) {
                knownContamination /=100;
            }
            std::cout << "..Contamination by normal cells set to:\t" << knownContamination*100<< "%\n";
    }
    if (knownContamination<0 && contaminationAdjustment == true) {
            std::cerr << "..contamination by normal cells should be a positive value\n";
            exit(0);
    }
    if (knownContamination==0 && contaminationAdjustment == true) {
        std::cout << "..FREEC is going to evaluate contamination by normal cells\n";
    }

    string pathToSamtools = (std::string)cf.Value("general","samtools","samtools");

	bool has_window = cf.hasValue("general","window");
    int window = (int)cf.Value("general","window",NA);


    bool has_coefficientOfVariation = cf.hasValue("general","coefficientOfVariation");
    float coefficientOfVariation = (float)cf.Value("general","coefficientOfVariation", 0.05);
    if (has_coefficientOfVariation && !has_window) {
            cout << "..Coefficient Of Variation set equal to "<< coefficientOfVariation<< "\n..it will be used to evaluate window size\n";
            if (coefficientOfVariation<=0) {
                cerr << "Error: 'coefficientOfVariation' must be positive\n";
		        cout << "..Since coefficientOfVariation' must be positive, FREEC will continue running with coefficientOfVariation=0.05\n";
                coefficientOfVariation=0.05;
            }
    } else if (has_coefficientOfVariation && has_window) {
            cout << "..Note, the Coefficient Of Variation won't be used since \"window\" = "<< window << " was set\n";
    } else if (!has_coefficientOfVariation && has_window) {
            cout << "..Window = "<< window << " was set\n";
    } else {
        cerr << "Error: 'coefficientOfVariation' or 'window' must be provided\n";
        cout << "..FREEC will use the coefficientOfVariation=0.05 to evaluate window size\n";
        coefficientOfVariation=0.05;
        has_coefficientOfVariation=true;
    }

    int step = (int)cf.Value("general","step", NA);
    if (step>0 && has_window && step < window) {
		cout << "..Step:\t" << step<< "\n";
    } else if (has_window) {
        step=window;
    } else if (step>0 && !has_window) {
        cerr << "Cannot set 'step' without 'window'\n";
        cout << "..Will ignore the value of step since window size is not provided\n";
        step = NA;
    }

	string outputDir = (std::string)cf.Value("general","outputDir",".");
	if ( access( outputDir.c_str(), 0 ) == 0 )    {
            struct stat status;
            stat( outputDir.c_str(), &status );

            if ( status.st_mode & S_IFDIR )   {
                cout << "..Output directory:\t" << outputDir << "\n";
            }
            else      {
                cerr << "Error: The path you entered for 'outputDir': "<< outputDir <<" is a file. It shoud be a directory" << endl;
                exit(-1);
            }
    }  else   {
            cerr << "Error: Path "<<outputDir<<" doesn't exist." << endl;
            exit(-1);
    }

	bool has_dirWithFastaSeq = cf.hasValue("general","chrFiles");
    string dirWithFastaSeq = (std::string)cf.Value("general","chrFiles","");
    if (has_dirWithFastaSeq) {
            if ( access( dirWithFastaSeq.c_str(), 0 ) == 0 )      {
                struct stat status;
                stat( dirWithFastaSeq.c_str(), &status );

                if ( status.st_mode & S_IFDIR )    {
                    cout << "..Directory with files containing chromosome sequences:\t" << dirWithFastaSeq << "\n";
                }
                else   {
                    cerr << "Error: The path you entered for 'dirWithFastaSeq': "<< dirWithFastaSeq <<" is a file. It shoud be a directory" << endl;
                    exit(-1);
                }
            }   else    {
                cerr << "Error: Path "<<dirWithFastaSeq<<" doesn't exist. Comment the line with 'chrFiles' if you use a precalculated GC-content profile or a control sample. Otherwise, set the correct path" << endl;
                exit(-1);
            }
    }

	int minimalTotalLetterCountPerPosition = round(float(cf.Value("BAF","minimalCoveragePerPosition", 0)));
	if (minimalTotalLetterCountPerPosition>0) {
        cout << "..will use a threshold of "<< minimalTotalLetterCountPerPosition <<" read(s) per SNP position to calculate beta allel frequency (BAF) values\n";
	}

	int minimalQualityPerPosition = int(cf.Value("BAF","minimalQualityPerPosition",0));
    int shiftInQuality =(int)cf.Value("BAF","shiftInQuality",0);
    if (minimalQualityPerPosition>0) {
        cout << "..will use a quality threshold of "<< minimalQualityPerPosition <<" to select nucleotides used in calculation of beta allel frequency (BAF) values\n";
        cout << "..will shift qualities by "<< shiftInQuality <<" when selecting nucleotides used in calculation of beta allel frequency (BAF) values\n";
        cout << "..Note, use shiftInQuality=33 for Sanger or Illumina 1.8+ format; shiftInQuality=64 for Illumina 1.3+\n";
        minimalQualityPerPosition += shiftInQuality;
    }

	bool has_sample_MateFile = cf.hasValue("sample","mateFile");
	bool has_sample_mateCopyNumberFile = cf.hasValue("sample","mateCopyNumberFile");
	std::string sample_MateFile = "";
	std::string sample_MateCopyNumberFile = "";
    string sample_inputFormat = std::string(cf.Value("sample","inputFormat",""));
    std::string sample_mateOrientation = (std::string)cf.Value("sample","mateOrientation","0");

    if (has_sample_MateFile) {
        sample_MateFile = std::string(cf.Value("sample","mateFile")) ;
		cout << "..Sample file:\t" << sample_MateFile << "\n";
        if (sample_inputFormat.compare("")==0) {
            cerr << "Error: You need to set the inputFormat to be avaible to read "<< sample_MateFile << "\n";
            cerr << "Available formats:SAM, BAM, pileup, Eland, BED, SOAP, arachne, psl (BLAT) and Bowtie\n";
            cerr << "FREEC works exclusively with 'inputFormat=pileup' when the user uses option [BAF]\n";
            exit (0);
        } else {
        		cout << "..Sample input format:\t" << sample_inputFormat << "\n";
        }
		if (sample_inputFormat.compare("BAM")==0 || sample_inputFormat.compare("bam")==0 || sample_inputFormat.compare("Bam")==0) {
            cout << "..will use this instance of samtools: '"<< pathToSamtools<<"' to read BAM files\n";
		}
    }
    if (has_sample_mateCopyNumberFile){
        sample_MateCopyNumberFile = std::string(cf.Value("sample","mateCopyNumberFile"));
        cout << "..Sample file with precalculated copy numbers:\t" << sample_MateCopyNumberFile << "\n";
    }
    if (!has_sample_MateFile && !has_sample_mateCopyNumberFile) {
         cerr << "Error: either \"mateFile\" or \"mateCopyNumberFile\" must be specified\n\n";
         exit(0);
    }

 	string myName = "";
    if (has_sample_MateFile) {
        myName = sample_MateFile;
    } else {
        myName = sample_MateCopyNumberFile;
    }

    bool has_control_MateFile = cf.hasValue("control","mateFile");
	bool has_control_mateCopyNumberFile = cf.hasValue("control","mateCopyNumberFile");
	std::string control_MateFile = "";
	std::string control_MateCopyNumberFile = "";
    string control_inputFormat = std::string(cf.Value("control","inputFormat",""));
    std::string control_mateOrientation = (std::string)cf.Value("control","mateOrientation","0");

    if (has_control_MateFile) {
        control_MateFile = std::string(cf.Value("control","mateFile")) ;
		cout << "..Control file:\t" << control_MateFile << "\n";
        if (control_inputFormat.compare("")==0) {
            cerr << "Error: You need to set the inputFormat to be avaible to read "<< control_MateFile << "\n";
            cerr << "Available formats:SAM, BAM, pileup, Eland, BED, SOAP, arachne, psl (BLAT) and Bowtie\n";
            cerr << "FREEC works exclusively with 'inputFormat=pileup' when the user uses option [BAF]\n";
            exit (0);
        } else {
        		cout << "..Input format for the control file:\t" << control_inputFormat << "\n";
        }
    }
    if (has_control_mateCopyNumberFile){
        control_MateCopyNumberFile = std::string(cf.Value("control","mateCopyNumberFile"));
        cout << "..Control file with precalculated copy numbers:\t" << control_MateCopyNumberFile << "\n";
    }

	bool isControlIsPresent = has_control_MateFile || has_control_mateCopyNumberFile;

	string controlName = "";
    if (has_control_MateFile) {
        controlName = control_MateFile;
    } else if (has_control_mateCopyNumberFile){
        controlName = control_MateCopyNumberFile;
    }

	bool sample_copyNumber_pileup_read = false;
	bool control_copyNumber_pileup_read = false;
	bool is_sample_pileup = (sample_inputFormat.compare("pileup") == 0 || sample_inputFormat.compare("SAMtools pileup") == 0);
	bool is_control_pileup = (control_inputFormat.compare("pileup") == 0 || control_inputFormat.compare("SAMtools pileup") == 0);

	bool has_BAF = cf.hasValue("BAF","SNPfile");
	if (has_BAF && !has_sample_MateFile) {
        cerr << "ERROR: you need to provide a 'mateFile' for the [sample] (in SAMtools pileup format) to be able to calculate BAF profiles with options [BAF] \n";
        exit (0);
	}

    if (has_BAF && !has_control_MateFile && isControlIsPresent) {
        cerr << "ERROR: you need to provide a 'mateFile' for the [control] (in SAMtools pileup format) to be able to calculate BAF profiles with options [BAF] and detect somatic CNAs and LOH\n";
        cerr << "..Otherwise, you may not to use the control data at all. Just comment or delete 'mateCopyNumberFile' in the [control] group of parameters\n";
        exit (0);
	}

	if (!is_sample_pileup && has_BAF) {
        cerr << "Error: to calculate BAF values, you need to provide mateFile in SAMtools pileup format\n";
        cout << "..since you mateFile is not in SAMtools pileup format, the BAF values will not be calculated\n";
        has_BAF=false;
	}
    string SNPinfoFile = std::string(cf.Value("BAF","SNPfile",""));

    bool ifTargeted = cf.hasValue("target","captureRegions");
    std::string targetBed = std::string(cf.Value("target","captureRegions",""));
    bool logLogNorm = false;
    //if (ifTargeted)logLogNorm=true;
    if (ifTargeted && !isControlIsPresent) {
        cerr << "ERROR: Currently you need to provide a control sample ('mateFile' or 'mateCopyNumberFile') when you analyze targeted sequencing data to eliminate capture bias. The GC-content bias is not the only bias in targeted sequencing\n";
        exit(0);
    }
    if (!has_window && ifTargeted) {
        cerr << "ERROR: At this point you need to profide window size when you work with a targeted sequencing experiment: option 'window' in group of parameters [general] in your config file\n";
        cerr << "choose a window for have at least several hundred reads per window\n";
        exit(0);
    }
    logLogNorm =bool(cf.Value("general","logLogNorm",logLogNorm)) ;
    //if (ifTargeted && !logLogNorm) {
    //    cerr << "Warning: I would recommend using logLogNorm=TRUE when working with targeted sequencing data\n";
    //}

    float minExpectedGC = float(cf.Value("general","minExpectedGC",0.35));
    float maxExpectedGC = float(cf.Value("general","maxExpectedGC",0.55));

	bool has_GCprofile = cf.hasValue("general","GCcontentProfile");
    std::string GCprofileFile = std::string(cf.Value("general","GCcontentProfile", ""));
    int forceGC = int(cf.Value("general","forceGCcontentNormalization",0));
    int intercept;
    bool isUseGC = false;
    if (!isControlIsPresent || has_BAF || forceGC) {
            if (!has_dirWithFastaSeq && !has_GCprofile) {
                cerr << "Error: with the current options, either 'chrFiles' or 'GCcontentProfile' must be set\n";
                exit(0);
            }
            isUseGC = true;
            if (ifTargeted) {
                if (forceGC==0) {
                    isUseGC = false;
                    cout <<"..Since you use targeted sequencing data, FREEC will use only control read counts to normalize copy number profiles.\n";
                    cout <<"....Set forceGCcontentNormalization=1 if you want to use GC-content normalization prior to control density normalization for targeted sequencing.\n";
                    cout <<"....However, with targeted sequencing, I would not recommend to use this option (forceGCcontentNormalization=1 or 2) since capture bias can be much stronger than GC-content bias\n";
                } else {
                    cout << "Warning: with targeted sequencing, I would not recommend to use forceGCcontentNormalization=1 or 2 since capture bias can be much stronger than GC-content bias\n";
                    cout <<"..I recommend you to set forceGCcontentNormalization=0 or comment this line in the config file\n";
                    cout << "..Continue anyway :-/\n";
                }
            }
    }

    if (isUseGC) {
        cout << "..minimal expected GC-content (general parameter \"minExpectedGC\") was set to "<< minExpectedGC<<"\n";
        cout << "..maximal expected GC-content (general parameter \"maxExpectedGC\") was set to "<< maxExpectedGC <<"\n";
        intercept = (int)(cf.Value("general","intercept", 1));
        if (intercept!=1) {
                cout << "Warning: I would advise using 'intercept=1' with your parameters\n";
        }
    }else {
        intercept = (int)(cf.Value("general","intercept", 0));
        if (intercept!=0) {
                cout << "Warning: I would advise using 'intercept=0' with your parameters\n";
        }
    }
    if (!ifTargeted && logLogNorm && isUseGC) {
        cerr << "Warning: will not use loglog-normalization since GC-content will be used\n";
        logLogNorm=false;
    }

    if (forceGC==2) { //will  Use GC, with intercept=0
         intercept = (int)(cf.Value("general","intercept", 0));
         if (intercept!=0) {
                cout << "Warning: I would advise using 'intercept=0' with your parameters\n";
         }
    }

    if(!cf.hasValue("general","chrLenFile")) {
        cerr <<"ERROR: you need to provide a file with chromosome lengths\n";
        exit(0);
    }
    std::string chrLenFile = (std::string)cf.Value("general","chrLenFile","");
    cout << "..File with chromosome lengths:\t" << chrLenFile << "\n";

    bool isMinMappabilitySet = cf.hasValue("general","minMappabilityPerWindow");
    minMappabilityPerWindow = double(cf.Value("general","minMappabilityPerWindow",0.85));
    if (isMinMappabilitySet && isUseGC) {
        cout << "..Using the minimal mappability of: "<< minMappabilityPerWindow <<"\n";
    } else if (isUseGC) {
        cout << "..Using the default minimal mappability value of "<< minMappabilityPerWindow <<"\n";
    }
    if (ifTargeted && !isUseGC) {
        cout << "..Mappability and GC-content won't be used\n";
        minMappabilityPerWindow = 0;
        cout << "..Control-FREEC won't use minimal mappability. All windows overlaping capture regions will be considered\n";
    }

    bool has_MapFile = cf.hasValue("general","gemMappabilityFile") ;
    string gemMapFile = std::string(cf.Value("general","gemMappabilityFile",""));

    if (cf.hasValue("general","uniqueMatch") && isUseGC) {
            string uMatch = std::string(cf.Value("general","uniqueMatch"));
            if (uMatch.compare("1")==0 || uMatch.compare("TRUE")==0 ||uMatch.compare("true")==0 ||uMatch.compare("True")==0) {
                if (!has_MapFile) {
                    cout << "Warning: FREEC set 'uniqueMatch=FALSE' since you did not provide a GEM mappability file ('gemMappabilityFile')";
                } else {
                    cout << "..Parameter uniqueMatch was set TRUE, will use "<< gemMapFile <<" for mappability information\n";
                    uniqueMatch = true;
                }
            }
    } else if (cf.hasValue("general","uniqueMatch")) {
        cout << "Warning: FREEC will not use option 'uniqueMatch' since FREEC is not going to use mappability or GC-content for normalization of copy number profiles\n";
    }
    if (!uniqueMatch) {
        cout << "..uniqueMatch = FALSE\n";
    }

    int ploidy;
    if(cf.hasValue("general","ploidy")){
        ploidy =round(float(cf.Value("general","ploidy")));
        cout << "..average ploidy set to "<<ploidy<<"\n";
    } else {
        ploidy = 2;
        cerr << "Warning: at this moment FREEC needs to know the genome average ploidy.. FREEC will set ploidy to "<<ploidy<<" (diploid genome)\n";
        cout << "Warning: at this moment FREEC needs to know the genome average ploidy.. FREEC will set ploidy to "<<ploidy<<" (diploid genome)\n";
		cout << "..Continue anyway\n";
    }

    int breakPointType= int(cf.Value("general","breakPointType",NORMALLEVEL));
    cout << "..break-point type set to "<<breakPointType<<"\n";

    bool noisyData = (bool)cf.Value("general","noisyData", "false");
    if ((!noisyData) && ifTargeted && has_BAF) {
        cout << "Warning: consider using '[general] noisyData=true' if you expect to have highly nonuniform coverage along the genome\n";
    } else if (noisyData && !has_BAF){
        cout << "Warning: Parameter '[general] noisyData=true' will not have effect since FREEC won't use BAF information to correct predicted copy numbers\n";
    }else if (noisyData &&  !ifTargeted ){
        cout << "Warning: I would not recommend using '[general] noisyData=true' for whole genome data; you can miss some real CNAs in this case\n";
    } else {
        cout << "..noisyData set to "<< noisyData<< "\n";
    }

    //if print -1 in the ratio files
    bool printNA = (bool)cf.Value("general","printNA", "true");

    int RCThresh = (int)cf.Value("general","readCountThreshold", 10);

    if(isControlIsPresent) {
        cout << "..minimal number of reads per window in the control sample is set to "<< RCThresh<< "\n";
    }

// createNames:

	char rsymb = '/';
	std::vector<std::string> elems = split(myName, '\\');
	if (elems.size()>1)
		rsymb = '\\';
	myName = elems.back();
	elems.clear();
	elems = split(myName, '/' );
	if (elems.size()>1)
		rsymb = '/';
	myName = elems.back();
	if (outputDir.compare("")!=0) {
		char lsymb = outputDir.at(outputDir.length()-1);
		if (lsymb != '/' && lsymb != '\\' ) {
			outputDir = outputDir+rsymb;
		}
			myName = outputDir+myName;
	}


	elems.clear();
	if (controlName.compare("")!=0) {
		elems = split(controlName, '\\');
		controlName = elems.back();
		elems.clear();
		elems = split(controlName, '/' );
		controlName = elems.back();
		controlName = outputDir+controlName;
		elems.clear();
	}

   //MAKE MINI PILEUP
    if (makeminipileup == true)
        {
        cout << "..Creating a Pileup file with SNP's in capture regions : \n";
        Minipileup minipileup;
        minipileup.makepileup(targetBed, sample_MateFile, sample_inputFormat, sample_mateOrientation, pathToSamtools);
        cout << "...-> Done!\n";
        }

    int given_ploidy = ploidy;
    std::vector <int> ploidies;
    std::vector <float> rmsscoresT;
    std::vector <float> rmsscoresC;
    //TRY OTHER PLOIDY
    if (tryOtherPloidy == true)
    {
        ploidies = vector<int>(4);
        for (int i = 0; i < 4; i++)
            {
            ploidies[i]=(ploidy-1+i);
            }
    }
    else
    {
    ploidies = vector<int>(1);
    ploidies[0] = ploidy;
    }

rmsscoresT = vector <float>(ploidies.size());
rmsscoresC = vector <float>(ploidies.size());
int test_ploidy = 0;
while (test_ploidy < ploidies.size())
{
    ploidy = ploidies[test_ploidy];
    cerr << "Running Control-FREEC with ploidy set to " << ploidy << "! \n";
    //READ SAMPLE DATA:

	GenomeCopyNumber sampleCopyNumber;
	sampleCopyNumber.setSamtools(pathToSamtools);

	GenomeCopyNumber controlCopyNumber;
	controlCopyNumber.setSamtools(pathToSamtools);

	SNPinGenome snpingenome;
	SNPinGenome snpingenomeControl;

	ThreadPoolManager* thrPoolManager = ThreadPoolManager::getInstance();
	ThreadPool* thrPool = NULL;
	if (has_BAF) {  //read the pileup files only once
	cerr << "Ca rentre! \n";
	  GenomeCopyNumberReadMateFileArgWrapper* readMateFileArg;

	  cout << "..will use SNP positions from "<< SNPinfoFile << " to calculate BAF profiles\n";
	  thrPool = thrPoolManager->newThreadPool("GenomeCopyNumber_readMateFile");

	  snpingenome.readSNPs(SNPinfoFile);

	  if (is_sample_pileup && !has_sample_mateCopyNumberFile && has_window) {
		std::cout << "avoid double pileup read: reading sample matefile\n";
		readMateFileArg = new GenomeCopyNumberReadMateFileArgWrapper(snpingenome, sample_MateFile, sample_inputFormat, minimalTotalLetterCountPerPosition, minimalQualityPerPosition, sampleCopyNumber, chrLenFile, window, step);
		thrPool->addThread(GenomeCopyNumber_readMateFile_wrapper, readMateFileArg);
		sample_copyNumber_pileup_read = true;
	  } else {
		readMateFileArg = new GenomeCopyNumberReadMateFileArgWrapper(snpingenome, sample_MateFile, sample_inputFormat, minimalTotalLetterCountPerPosition, minimalQualityPerPosition);
		thrPool->addThread(GenomeCopyNumber_readMateFile_wrapper, readMateFileArg);
		sample_copyNumber_pileup_read = false;
	  }

	  if (isControlIsPresent) {
		snpingenomeControl.setSNPChr(snpingenome.getSNPChr());

		if (is_control_pileup && !has_control_mateCopyNumberFile && has_window) {
		  std::cout << "avoid double pileup read: reading control matefile\n";
		  readMateFileArg = new GenomeCopyNumberReadMateFileArgWrapper(snpingenomeControl, control_MateFile, control_inputFormat, minimalTotalLetterCountPerPosition, minimalQualityPerPosition, controlCopyNumber, chrLenFile, window, step);
		  thrPool->addThread(GenomeCopyNumber_readMateFile_wrapper, readMateFileArg);
		  control_copyNumber_pileup_read = true;
		} else {
		  readMateFileArg = new GenomeCopyNumberReadMateFileArgWrapper(snpingenomeControl, control_MateFile, control_inputFormat, minimalTotalLetterCountPerPosition, minimalQualityPerPosition);
		  thrPool->addThread(GenomeCopyNumber_readMateFile_wrapper, readMateFileArg);
		  control_copyNumber_pileup_read = false;
		}
	  }
      thrPool->run();
      delete thrPool;

	}


    if((ifTargeted) && (window == 0) && (step == 0))
       {
       cout << "FREEC will work on exome sequencing data! \n";
       if (has_sample_mateCopyNumberFile) {
            sampleCopyNumber.readCopyNumber(sample_MateCopyNumberFile);
       } else {
            if (!sample_copyNumber_pileup_read && has_window) {
                    sampleCopyNumber.readCopyNumber(sample_MateFile, sample_inputFormat, sample_mateOrientation,chrLenFile, window, step, targetBed, ifTargeted);
            } else if (!sample_copyNumber_pileup_read && !has_window) {
                sampleCopyNumber.readCopyNumber(sample_MateFile, sample_inputFormat, sample_mateOrientation,chrLenFile, coefficientOfVariation);
            }
            sampleCopyNumber.printCopyNumber(myName+"_sample.cpn");
        }
       }
       else {
        if (step != NA) {
            sampleCopyNumber.setStep(step);
        }

        if (has_sample_mateCopyNumberFile) {
            sampleCopyNumber.readCopyNumber(sample_MateCopyNumberFile);
            step = sampleCopyNumber.getStep();
        }  else {
            if (!sample_copyNumber_pileup_read && has_window) {
                    sampleCopyNumber.readCopyNumber(sample_MateFile, sample_inputFormat, sample_mateOrientation,chrLenFile, window, step);
            } else if (!sample_copyNumber_pileup_read && !has_window) {
                sampleCopyNumber.readCopyNumber(sample_MateFile, sample_inputFormat, sample_mateOrientation,chrLenFile, coefficientOfVariation);
                step = sampleCopyNumber.getWindowSize(); //in this case step=windowSize
            }
            sampleCopyNumber.printCopyNumber(myName+"_sample.cpn");
        }
        window = sampleCopyNumber.getWindowSize();
        cout << "..Window size:\t"<< window << "\n";
        if (step == NA) {
            step= window;
        }
        has_window = true; //now we know window size and even step!
        sampleCopyNumber.setSex(sex);
        }

    //READ CONTROL DATA:


if (isControlIsPresent) {
	if (has_control_mateCopyNumberFile) {
		controlCopyNumber.readCopyNumber(control_MateCopyNumberFile);
	} else {
        if (!control_copyNumber_pileup_read ) {
				  controlCopyNumber.readCopyNumber(control_MateFile, control_inputFormat, control_mateOrientation, chrLenFile, window, step, targetBed, ifTargeted);
		}
      controlCopyNumber.printCopyNumber(controlName+"_control.cpn");
	}
   controlCopyNumber.setSex(sex);
}

    //if it is a TARGETED resequencing experiment, delete all info outside of the target regions
	if(ifTargeted && (step != 0 || window != 0)) {
        cout << "..FREEC will take into account only regions from "<<targetBed<<"\n";
        int minRegion = sampleCopyNumber.focusOnCapture(targetBed);
        if (teloCentroFlanks>minRegion) {
            teloCentroFlanks = minRegion;
            cout << "..telocenromeric set to "<<teloCentroFlanks<<" since it is the minimal length of capture regions\n";
        }
        controlCopyNumber.focusOnCapture(targetBed);
	}

    //READ GC-CONTENT:

	if (isUseGC) {//then read CG-content.
		cout << "..using GC-content to normalize copy number profiles\n";
        if (has_GCprofile) {			// a file with CG-content already exists
			int stepGC = sampleCopyNumber.readCGprofile(GCprofileFile);
			if (step!=stepGC) {
			    cerr << "Error: Uncorrect window size in the GC-content profile. FREEC will need to recalculate it. You must provide a path to chromosome files, option \"chrFiles\"\n";
                exit (0);
			}
        } else {// has_dirWithFastaSeq is true
			sampleCopyNumber.fillCGprofile(dirWithFastaSeq);
			GCprofileFile = outputDir+"GC_profile.cnp";
			if (!has_MapFile)
                sampleCopyNumber.printCGprofile(GCprofileFile); //if has_MapFile will print out GC-content later
        }
        if (has_MapFile) {  //read mappability file
            sampleCopyNumber.readGemMappabilityFile(gemMapFile);
            //rewrite GC-profile with mappability as the last (5th) colomn
            GCprofileFile = outputDir+"GC_profile.cnp";
            sampleCopyNumber.printCGprofile(GCprofileFile);
            cout << "..Mappability track from "<< gemMapFile <<" has been added to "<< GCprofileFile <<"\n";
        }
	}

	if (isControlIsPresent && isUseGC) {//then read CG-content and associate it with the control data.
		cout << "..using GC-content to normalize the control profile\n";
		controlCopyNumber.readCGprofile(GCprofileFile); //the file with CG-content already exists
        if (ifTargeted && window != 0) {
                controlCopyNumber.focusOnCapture(targetBed); // to mask again averything which is not in the capture
                sampleCopyNumber.focusOnCapture(targetBed); //TODO: Check whether it is needed. can get here only if "forceGCwhenControlIsPresent>0"
        }
	}
    //remove window with read count less than RCThresh from the analysis:
    if (isControlIsPresent)
        sampleCopyNumber.removeLowReadCountWindows(controlCopyNumber,RCThresh);

    //NORMALIZE READ COUNT:
    sampleCopyNumber.setPloidy(ploidy);
	if (isControlIsPresent) {
		//check if window size is the same for the Control and Sample
		if (sampleCopyNumber.getWindowSize() != controlCopyNumber.getWindowSize()) {
			cerr << "\nError: the window length is different for sample and control data\n\tPlease check parameters and input files!\n\n";
			return -1;
		}
		if ((!forceGC && !has_BAF) || (ifTargeted&&forceGC!=1)) { //normalize sample density with control density
            sampleCopyNumber.calculateRatio(controlCopyNumber, degree,intercept,logLogNorm);
		} else { //forceGC != 0
		    if (forceGC==1) { //normalize first Sample and Control, and then calculate the ratio
                controlCopyNumber.setPloidy(2); // normal genome has ploidy=2!!!
                cout << "..Set ploidy for the control genome equal to "<< 2 << "\n";
                sampleCopyNumber.calculateRatioUsingCG(degree, intercept,minExpectedGC,maxExpectedGC);
                //exit(0);
                controlCopyNumber.calculateRatioUsingCG(degree, intercept,minExpectedGC,maxExpectedGC);
                sampleCopyNumber.calculateRatioUsingCG(controlCopyNumber);
            } else if (forceGC==2) {  //calculate the ratio , normalize for GC
                sampleCopyNumber.calculateRatio(controlCopyNumber, degree,intercept,logLogNorm);
                sampleCopyNumber.recalculateRatioUsingCG(8, 1,minExpectedGC,maxExpectedGC); //try higher values of polynomial's degree
		    }
		}
        if(has_BAF && forceGC!=1 && !ifTargeted) { //calculateRatioUsingCG
            if (intercept != 1) cout << "Warning: Again, I would advise using 'intercept = 1' with your parameters\n";

            if (forceGC==0) { //otherwise, already calculated for the Sample
                sampleCopyNumber.calculateRatioUsingCG(degree, intercept,minExpectedGC,maxExpectedGC);
		    }
            controlCopyNumber.setPloidy(2); //normal genome has ploidy=2
            cout << "..Set ploidy for the control genome equal to "<< 2 << "\n";
            controlCopyNumber.calculateRatioUsingCG(degree, intercept,minExpectedGC,maxExpectedGC);

        }
        if (ifTargeted && has_BAF && forceGC!=1) {
            cout << "Warning: Control-FREEC will assume that there is not gains and losses in the target regions in the control genome\n";
            controlCopyNumber.setPloidy(2); //normal genome has ploidy=2
            cout << "..Set copy number in the control genome equal to "<< 2 << "\n";
            controlCopyNumber.setAllNormal();
        }

	} else { // no Control present
        sampleCopyNumber.setPloidy(ploidy);
		sampleCopyNumber.calculateRatioUsingCG(degree, intercept,minExpectedGC,maxExpectedGC);
	}
	cout << "..Copy number profile normalization -> done\n";
	//segmentation:

    if (knownContamination>0) {
        cout << "..Recalculating copy number profiles using known value of contamination by normal cells:\n";
        cout << ".."<< knownContamination*100 <<"%\n";
        sampleCopyNumber.recalculateRatio(knownContamination);
        sampleCopyNumber.setNormalContamination(knownContamination);
    }

	thrPool = thrPoolManager->newThreadPool("GenomeCopyNumber_calculateBreakpoint");
	GenomeCopyNumberCalculateBreakpointArgWrapper* bkpArg = new GenomeCopyNumberCalculateBreakpointArgWrapper(sampleCopyNumber, breakPointThreshold, breakPointType);
	thrPool->addThread(GenomeCopyNumber_calculateBreakpoint_wrapper, bkpArg);
	if (has_BAF && isControlIsPresent) {
	  bkpArg = new GenomeCopyNumberCalculateBreakpointArgWrapper(controlCopyNumber, breakPointThreshold, breakPointType);
	  thrPool->addThread(GenomeCopyNumber_calculateBreakpoint_wrapper, bkpArg);
	}
	thrPool->run();
	delete thrPool;

	cout << "..calculate median values\n";
	cout << std::flush;
	sampleCopyNumber.calculateCopyNumberMedians(minCNAlength,0);
	//sampleCopyNumber.calculatePloidy(minCNAlength);
	sampleCopyNumber.recalcFlanks(teloCentroFlanks, 3);
	//sampleCopyNumber.deleteFlanks(TELO_CENTRO_FLANCS);
	cout << "..annotate copy numbers\n";
	cout << std::flush;
	sampleCopyNumber.calculateCopyNumberProbs_and_genomeLength(breakPointType);
	if (contaminationAdjustment && knownContamination==0) {
	    cout <<"..Evaluating possible contamination..\n";
		cout << std::flush;
		float contamValue=sampleCopyNumber.evaluateContamination();
		cout << "..Identified contamination by normal cells: " << contamValue*100 << "%\n";
		cout << std::flush;
		if (contamValue>0) {
		    cout << "..Recalculating copy number profiles..\n";
			cout << std::flush;
            sampleCopyNumber.recalculateRatio(contamValue);
            knownContamination = contamValue;
            sampleCopyNumber.setNormalContamination(knownContamination);
            cout << "..Recalculate breakpoints\n";
			cout << std::flush;
            sampleCopyNumber.calculateBreakpoints(breakPointThreshold,breakPointType);
            cout << "..Recalculate median values\n";
            sampleCopyNumber.calculateCopyNumberMedians(minCNAlength,0);
			cout << std::flush;
            //sampleCopyNumber.deleteFlanks(TELO_CENTRO_FLANCS);
            sampleCopyNumber.recalcFlanks(teloCentroFlanks, 3);
            cout << "..Reannotate copy numbers\n";
            sampleCopyNumber.calculateCopyNumberProbs_and_genomeLength(breakPointType);
			cout << std::flush;
            //sampleCopyNumber.printRatio(myName+"_ratio.txt");
		}
	}
	if (has_BAF) {
        breakPointThreshold = 0.8;
        if (ifTargeted)
            breakPointThreshold = 1.6;

		thrPool = thrPoolManager->newThreadPool("SNPinGenome_perform");

		SNPinGenomePerformArgWrapper* snpArg = new SNPinGenomePerformArgWrapper(snpingenome, sample_MateFile, sample_inputFormat, minimalTotalLetterCountPerPosition,minimalQualityPerPosition, noisyData,sampleCopyNumber, breakPointThreshold, breakPointType, minCNAlength, "Sample");
		thrPool->addThread(SNPinGenome_perform_wrapper, snpArg);

        //the same for the control sample:
        if (isControlIsPresent) {
			snpArg = new SNPinGenomePerformArgWrapper(snpingenomeControl, control_MateFile, control_inputFormat, minimalTotalLetterCountPerPosition,minimalQualityPerPosition, noisyData, controlCopyNumber, breakPointThreshold, breakPointType, minCNAlength, "Control");
			thrPool->addThread(SNPinGenome_perform_wrapper, snpArg);
		}

		thrPool->run();
		delete thrPool;

	    sampleCopyNumber.printBAF(myName+"_BAF.txt",snpingenome);

        if (isControlIsPresent) {
            sampleCopyNumber.calculateSomaticCNVs(controlCopyNumber.getCNVs(),controlCopyNumber.getPloidy());

            controlCopyNumber.printBAF(myName+"_normal_BAF.txt",snpingenomeControl);
            controlCopyNumber.printRatio(myName+"_normal_ratio.txt",0,printNA);

            if (ifBedGraphOutPut) {
                controlCopyNumber.printRatio(myName+"_normal_ratio.BedGraph",1,printNA);
            }

            controlCopyNumber.printCNVs(myName+"_normal_CNVs");
        }
    }

	sampleCopyNumber.printRatio(myName+"_ratio.txt",0,printNA);
	if (ifBedGraphOutPut) {
            sampleCopyNumber.printRatio(myName+"_ratio.BedGraph",1,printNA);
    }
	sampleCopyNumber.printCNVs(myName+"_CNVs");
    test_ploidy++;
}

cerr << "\n";
cerr << "Ploidy" << "\t" << "RMStumor" << "\n"; // << "RMS score Control" << "\n";
for (int  i=0; i < ploidies.size(); i++)
{
cerr << ploidies[i] << "\t" << rmsscoresT[i] << "\n"; // << rmsscoresC[i] << "\n";
}
	return 0;
}
