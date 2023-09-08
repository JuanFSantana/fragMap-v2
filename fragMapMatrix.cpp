#include <iostream>
#include <string>
#include <map>
#include <fstream> 
#include <nlohmann/json.hpp>


void establishEvents(std::ifstream& reads_bed_file, std::map<int, std::map<int, double>>& intervals, const std::string& analysis_type, const int& userStart, const int& userEnd, double& correctionFactor, int& userFragStart, int& userFragEnd);
void computeCounts(std::map<int, std::map<int, double>>& intervals, std::map<int, std::map<int, double>>& intervalsCounts);
void fillInMissingPositions(std::map<int, std::map<int, double>>& intervalsCounts, std::map<int, std::map<int, double>>& intervalsCountsFinal, const int& userStart, const int& userEnd, int& userFragStart, int& userFragEnd);
void outputTheFile(std::ofstream& outputFile, std::map<int, std::map<int, double>>& intervalsCountsFinal, const std::string& analysis_type, const std::string& outputName); 

int main(int argc, char* argv[]) {
    // make sure that 6 arguments are passed
    if (argc != 8) {
        std::cerr << "Usage: " << argv[0] << " <json data> <type of analysis(centers or full frag)> <output name> <start of regions coordinate (-500)> <end of regions coordinate (500)> <start of frag lengths (20)> <end of frag lengths (400)>" << std::endl << std::endl;
        exit(EXIT_FAILURE);
        return 1;
    }
    // check for valid type_of_overlap option and establish type of analysis; declaring bedFilePath to be opened, jsonData to be parsed 
    std::string receivedJsonString = argv[1];
    std::string analysis_type = argv[2];
    std::string outputName = argv[3];
    int userStart = std::stoi(argv[4]);
    int userEnd = std::stoi(argv[5]);
    int userFragStart = std::stoi(argv[6]);
    int userFragEnd = std::stoi(argv[7]);

    nlohmann::json data = nlohmann::json::parse(receivedJsonString);

    // create vector of Interval objects
    std::map<int, std::map<int, double>> intervals;
    try {
        for (auto& element : data.items()) {
            std::string bedFilePath = element.key();
            
            // Check if the value is a number
            if (element.value().is_number()) {
                double correctionFactor = element.value().get<double>();
                std::ifstream reads_bed_file(bedFilePath);
                if (!reads_bed_file.is_open()) {
                    std::cerr << "Unable to open file: " << bedFilePath << std::endl;
                    continue;
                }
                establishEvents(reads_bed_file, intervals, analysis_type, userStart, userEnd, correctionFactor, userFragStart, userFragEnd);
                reads_bed_file.close();

            } else {
                std::cerr << "Unexpected value for key '" << bedFilePath << "'. Expected a number." << std::endl;
            }
        }
        } catch (const nlohmann::json::exception& e) {
            std::cerr << "JSON error: " << e.what() << std::endl;
        } catch (const std::exception& e) {
            std::cerr << "Standard exception: " << e.what() << std::endl;
        } catch (...) {
            std::cerr << "Unknown error occurred!" << std::endl;
        }

    
    std::map<int, std::map<int, double>> intervalsCounts;
    computeCounts(intervals, intervalsCounts); 
    // done with intervals, so clear it
    intervals.clear();
    // make final map with all positions filled in
    std::map<int, std::map<int, double>> intervalsCountsFinal;  
    fillInMissingPositions(intervalsCounts, intervalsCountsFinal, userStart, userEnd, userFragStart, userFragEnd);
    // done with intervalsCounts, so clear it
    intervalsCounts.clear();
    // output file
    std::ofstream outputFile;
    outputTheFile(outputFile, intervalsCountsFinal, analysis_type, outputName);

    return 0;
}

void establishEvents(std::ifstream& reads_bed_file, std::map<int, std::map<int, double>>& intervals, const std::string& analysis_type, const int& userStart, const int& userEnd, double& correctionFactor, int& userFragStart, int& userFragEnd){
    // read each line of the reads bed file
    // declare variables to store bed file information from bedtools output
    std::string regionChr; // col1
    int regionStart; // col2
    int regionEnd; // col3
    std::string regionGene; // col4
    std::string additionalInfo1; // col5
    std::string regionStrand; // col6      
    std::string readsChr; // col7
    int readsStart; // col8
    int readsEnd; // col9
    std::string readsID; // col 10 
    int readsQual; // col 11
    std::string readsStrand; // col12 

    while (reads_bed_file >> regionChr >> regionStart >> regionEnd >> regionGene >> additionalInfo1 >> regionStrand >> readsChr >> readsStart >> readsEnd >> readsID >> readsQual >> readsStrand) {
        int readSize = (readsEnd - readsStart) + 1;
        if (readSize >= userFragStart && readSize <= userFragEnd){
            int readBasePositionStart;
            int readBasePositionEnd;
            int positionTSS = (regionStart + regionEnd) / 2;
            int newStart = -99999;
            int newEnd = -99999;
            // calculate read position relative to TSS in the middle of the region
            if (analysis_type == "centers"){
                if (regionStrand == "+"){
                    readBasePositionStart = ((readsStart + readsEnd) / 2) - positionTSS;  
                    readBasePositionEnd = readBasePositionStart;
                }

                else if (regionStrand == "-"){
                    readBasePositionStart = positionTSS - (readsStart + readsEnd) / 2; 
                    readBasePositionEnd = readBasePositionStart;
                }

                // if the postion is 0 or positive, add 1 to the position to reflect starting at +1
                if (readBasePositionStart >= 0){
                    readBasePositionStart += 1;
                }
                if (readBasePositionEnd >= 0){
                    readBasePositionEnd += 1;
                }

                if (readBasePositionStart >= userStart && readBasePositionStart <= userEnd && readBasePositionEnd >= userStart && readBasePositionEnd <= userEnd){
                    newStart = readBasePositionStart;
                    newEnd = readBasePositionEnd;
                } else{
                    continue;
                }
            }

            else if (analysis_type == "full"){
                if (regionStrand == "+"){
                    readBasePositionStart = readsStart - positionTSS;
                    readBasePositionEnd = readsEnd - positionTSS;
                }

                else if (regionStrand == "-"){
                    readBasePositionStart = positionTSS - readsEnd;
                    readBasePositionEnd = positionTSS - readsStart;
                }

                // if the postion is 0 or positive, add 1 to the position to reflect starting at +1
                if (readBasePositionStart >= 0){
                    readBasePositionStart += 1;
                }
                if (readBasePositionEnd >= 0){
                    readBasePositionEnd += 1;
                }
                newStart = std::max(readBasePositionStart, userStart);
                newEnd = std::min(readBasePositionEnd, userEnd);
            }

            // 4 types of reads to count:
            // 1. both left and right position of read are outside of user defined start and end but read overlaps with user defined start and end
            // 2. left position of read is less than start but right position of read is inside user defined start and end
            // 3. left position of read is inside user defined start and end but right position of read is greater than end
            // 4. both left and right position of read are inside user defined start and end
            // addto intervals
            if (newStart != -99999 && newEnd != -99999){
                if (intervals[readSize].find(newStart) == intervals[readSize].end()) {
                        intervals[readSize][newStart] = 1 * correctionFactor; // add 1 to the count and normalize by the correction factor

                } else {
                    intervals[readSize][newStart] += 1 * correctionFactor; 
                }

                int newEndNotZero; // for cases where newStart is -1, newEnd is 0, we want to add 1 to newEnd
                if ((newEnd+1) == 0){
                    newEndNotZero = 1;
                }
                else{
                    newEndNotZero = newEnd+1;
                }
                if (intervals[readSize].find(newEnd) == intervals[readSize].end()) {

                    intervals[readSize][newEndNotZero] = -1 * correctionFactor;
                } else {
                    intervals[readSize][newEndNotZero] -= 1 * correctionFactor;
                }
            }
        }
    }
}

void computeCounts(std::map<int, std::map<int, double>>& intervals, std::map<int, std::map<int, double>>& intervalsCounts) {
    // iterate through each key-value pair in the map

    for (const auto& [readLength, intervalPairs]: intervals){
        // initialize current count
        double currentCount = 0;
        // iterate through each pair in the vector
        for (const auto& [position, change]: intervalPairs){
            // update the current count
            currentCount += change;
            // update the result in intervalsCounts
            intervalsCounts[readLength][position] = currentCount;
        }
    }
}

void fillInMissingPositions(std::map<int, std::map<int, double>>& intervalsCounts, std::map<int, std::map<int, double>>& intervalsCountsFinal, const int& userStart, const int& userEnd, int& userFragStart, int& userFragEnd){
    // iterate through each key-value pair in the map
    for (const auto& [readLength, intervalPairs]: intervalsCounts) {
        // initialize vector with seen positions
        std::vector<int> seenPositions;
        // first loop: populate seenPositions
        for (const auto& [position, count]: intervalPairs) {
            seenPositions.push_back(position);
        }
        // second loop: iterate through seenPositions to fill missing positions
        for (size_t eachPosition = 0; eachPosition < seenPositions.size() - 1; eachPosition++) {
            for (int absentPositions = seenPositions[eachPosition] + 1; 
                 absentPositions < seenPositions[eachPosition + 1]; absentPositions++) {
                if (absentPositions != 0) {  // skip if absentPositions is 0
                    intervalsCountsFinal[readLength][absentPositions] = intervalsCounts[readLength][seenPositions[eachPosition]];
                }
            }
        }

        // copy existing positions to the final map
        for (const auto& [position, count] : intervalPairs) {
            if (position <= userEnd){ // skip last position if is greater than userEnd
                intervalsCountsFinal[readLength][position] = count;
            }
        }
    
    }
    // fill in any missing fragSizes and positions with 0
    for (int i = userFragStart; i <= userFragEnd; i++){
        auto& innerMap = intervalsCountsFinal[i];  // This will create the inner map if it doesn't exist
        for (int j = userStart; j <= userEnd; j++){
            if (j != 0 && innerMap.find(j) == innerMap.end()) {
                innerMap[j] = 0.0;
            }
        }
    }
}

void outputTheFile(std::ofstream& outputFile, std::map<int, std::map<int, double>>& intervalsCountsFinal, const std::string& analysis_type, const std::string& outputName){
    outputFile.open(outputName+".bed");
    // create lines
    std::string header = "\t";
    for (const auto& [readLength, intervalPairs]: intervalsCountsFinal){
        for (const auto& [position, count]: intervalPairs){
            header +=  std::to_string(position) + "\t";
        }
        // remove last tab
        header.pop_back();
        header = header + "\n";
        break;
    }

    outputFile << header;

    for (const auto& [readLength, intervalPairs]: intervalsCountsFinal){
        std::string size = std::to_string(readLength);
        std::string line = size+"\t";

        for (const auto& [position, count]: intervalPairs){
            line += std::to_string(count) + "\t";
        }
        // remove last tab
        line.pop_back();
        line = line + "\n";
        outputFile << line;
    }

    outputFile.close();
}
