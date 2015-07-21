#pragma once
#ifndef MINIPILEUP_H
#define MINIPILEUP_H

#include <iostream>
#include <stdlib.h>
#include <fstream>
#include <vector>
#include "myFunc.h"


class Minipileup
{
    public:
        Minipileup();
        virtual ~Minipileup();
        void makepileup(std::string targetBed, std::string const& mateFileName ,std::string const& inputFormat, std::string const& matesOrientation, std::string pathToSamtools);
        void intersectWithBedtools();
        float calculateFlankLength(std::string const& mateFileName ,std::string const& inputFormat, std::string const& matesOrientation, std::string pathToSamtools);
        void calculateNewBoundaries(std::string targetBed, int flanks);
        void printNewCaptureRegions();
        void createPileUpFile();
        std::vector <int> coordinates_;
        std::vector <int> ends_;
        int exons_Count;
        int exons_Counttmp;

    private:
        int length_;
        std::string chromosome_;
        std::vector <std::string> coordinatesTmp_;
        std::vector <std::string> endsTmp_;
        std::vector <std::string> chr_names;
        std::vector <int> length_with_flanks;
        std::vector <std::string> length_with_flanksTmp;
        std::vector <std::string> strand;
        std::vector <std::string> ref_name;

};

#endif // MINIPILEUP_H
