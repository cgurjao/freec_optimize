#include "Minipileup.h"

using namespace std;

Minipileup::Minipileup()
{
}

Minipileup::~Minipileup()
{
}

void Minipileup::makepileup(std::string targetBed, std::string const& mateFileName ,std::string const& inputFormat, std::string const& matesOrientation, std::string pathToSamtools)
{
    int flanks = calculateFlankLength(mateFileName, inputFormat, matesOrientation, pathToSamtools);
    calculateNewBoundaries(targetBed, flanks);
    printNewCaptureRegions();
    intersectWithBedtools();
    createPileUpFile();
    /*FILE *stream;
    string command = "/bioinfo/local/build/BEDTools_2.21.0/bin/bedtools intersect -a /data/tmp/cgurjao/Data_NB/1000G_omni2.5.hg19.sites.vcf -b /data/tmp/cgurjao/Data_NB/New_Capture_Regions_for_SNP.bed > /data/tmp/cgurjao/Results/Results_temp/SNP_in_exons_and_flanks.bed";
    stream = popen(command.c_str(), "w");
    pclose(stream);*/
}


float Minipileup::calculateFlankLength(std::string const& mateFileName, std::string const& inputFormat_str, std::string const& matesOrientation_str, std::string pathToSamtools_)
{
        float read_Size = 0;
        std::ifstream fileMates (mateFileName.c_str());
        vector<float> insertSizeVector;
        string line;
        if (!fileMates.is_open())
            {
            cerr << "Error: unable to open "+mateFileName+"\n" ;
            exit(-1);
            }

        #ifdef PROFILE_TRACE
            time_t t0 = time(NULL);
        #endif

            MateOrientation matesOrientation = getMateOrientation(matesOrientation_str);
            InputFormat inputFormat;
            char* line_buffer;
            FILE *stream;
            char buffer[MAX_BUFFER];
            string command;
            string myInputFormat=inputFormat_str;
            if (mateFileName.substr(mateFileName.size()-3,3).compare(".gz")!=0) {
                    command = pathToSamtools_ + " view "+mateFileName;
                    myInputFormat="sam";       //will try to use existing samtools
            }
            inputFormat = getInputFormat(myInputFormat);
            stream =
        #ifdef _WIN32
            _popen(command.c_str(), "r");
        #else
            popen(command.c_str(), "r");
        #endif
        int j = 0;
        while (((line_buffer = getLine(buffer, MAX_BUFFER, stream, line)) != NULL))
        {
        if (line_buffer[0] == '@')
            return 0;

        string chr1,chr2;
        char* strs[32];
        unsigned int strs_cnt = split((char*)line_buffer, '\t', strs);
        if (strs_cnt > 3) {
            string sequence = strs[9];
            read_Size += sequence.size();
            }
        j++;
        }
        float average_read_Size = read_Size/j;
        float flanks = average_read_Size/2;
        #ifdef _WIN32
				_pclose(stream);
		#else
				pclose(stream);
		#endif
        return flanks;
}

void Minipileup::calculateNewBoundaries(std::string targetBed, int flanks)
{
        std::string const& captureFile = targetBed ;
        ifstream file (captureFile.c_str());
        if (file.is_open())
            {
            std::string line;
            exons_Count = 0;
            exons_Counttmp = 0;
            while (std::getline(file,line))
                {
                    exons_Count++;
                }
            int l=0;
            file.clear();
            file.seekg(0);
            while (std::getline(file,line))
                {
                    int i = 0;
                    while (line[i] != 'r')
                        {
                        i++;
                        }
                        l++;
                }
            file.clear();
            file.seekg(0);
            l = 0;
            while (std::getline(file,line))
                {
                exons_Counttmp++;
                l++;
                }

            length_ = exons_Counttmp;
            chr_names = vector<string>(exons_Counttmp);
            coordinatesTmp_ = vector<string>(exons_Counttmp);
            endsTmp_ = vector<string>(exons_Counttmp);
            length_with_flanksTmp = vector<string>(exons_Counttmp);
            strand = vector<string>(exons_Counttmp);
            ref_name = vector<string>(exons_Counttmp);

            l=0;
            file.clear();
            file.seekg(0);

            while (std::getline(file,line) && l < exons_Counttmp)
                {
                    int i = 0;
                    while (line[i] != '\t')
                        {
                        chr_names[l] += line[i];
                        i++;
                        }
                        i++;
                    while (line[i] != '\t')
                        {
                        coordinatesTmp_[l] += line[i];
                        i++;
                        }
                        i++;
                    while (line[i] != '\t')
                        {
                        endsTmp_[l] += line[i];
                        i++;
                        }
                        i++;
                    while (line[i] != '\t')
                        {
                        ref_name[l] += line[i];
                        i++;
                        }
                        i++;
                    while (line[i] != '\t')
                        {
                        length_with_flanksTmp[l] += line[i];
                        i++;
                        }
                        i++;
                    strand[l] += line[i];
                    l++;
                }
            }
            else
            {
            cerr << "Error: Unable to open file "+captureFile+"\n";
            exit(-1);
            }

        coordinates_ = vector<int>(exons_Count);
        ends_ = vector<int>(exons_Count,0);
        length_with_flanks = vector<int>(exons_Count);
        for (int i = 0; i<exons_Counttmp; i++)
            {
            ends_[i] = atoi(endsTmp_[i].c_str()) + flanks ;
            }
        for (int i = 0; i<exons_Counttmp; i++)
            {
            coordinates_[i] = atoi(coordinatesTmp_[i].c_str()) - flanks;
            }
        for (int i = 0; i<exons_Counttmp; i++)
            {
            length_with_flanks[i] = atoi(length_with_flanksTmp[i].c_str()) + 2*flanks;
            }
}

void Minipileup::printNewCaptureRegions()
{
  ofstream myfile;
  myfile.open ("/data/tmp/cgurjao/Results/Results_temp/New_Capture_Regions_for_SNP.bed");
  for (int i = 0; i < exons_Count; i++)
  {
    myfile << chr_names[i] << "\t" << coordinates_[i] << "\t" << ends_[i] << "\t" << ref_name[i] << "\t" <<  length_with_flanks[i] << "\t" << strand[i] << "\n";
  }
  myfile.close();
}

void Minipileup::intersectWithBedtools()
{
    FILE *stream;
    string command = "/bioinfo/local/build/BEDTools_2.21.0/bin/bedtools intersect -a /data/tmp/cgurjao/Data_NB/1000G_omni2.5.hg19.sites.vcf -b /data/tmp/cgurjao/Results/Results_temp/New_Capture_Regions_for_SNP.bed > /data/tmp/cgurjao/Results/Results_temp/SNP_in_exons_and_flanks.bed";
    stream = popen(command.c_str(), "w");
    pclose(stream);
}

void Minipileup::createPileUpFile()
{
    /*FILE *stream;
    string command = "/bioinfo/local/build/BEDTools_2.21.0/bin/bedtools intersect -a /data/tmp/cgurjao/Data_NB/1000G_omni2.5.hg19.sites.vcf -b /data/tmp/cgurjao/Results/Results_temp/New_Capture_Regions_for_SNP.bed > /data/tmp/cgurjao/Results/Results_temp/SNP_in_exons_and_flanks.bed";
    stream = popen(command.c_str(), "w");
    pclose(stream);*/
}
