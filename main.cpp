#include <iostream>
#include <cstdlib>
#include "constants.h"
#include "linkage_group_DH.h"
#include "genetic_map_DH.h"
#include "genetic_map_RIL.h"

int main (int argc, char * const argv[]) {

    if (argc != 3) {
        cout << "mstmap v1.0 (2008, University of California, Riverside)" << endl;
        cout << "\t Usage: ./mstmap.exe input_file output_file" << endl;
	cout << "\t See example.txt for an example input file" << endl;
	cout << "\t Questions: email Stefano Lonardi <stelo@cs.ucr.edu>" << endl;
	cout << "\t Citation: Y.Wu, P.Bhat, T.J.Close, S.Lonardi" << endl;
	cout << "\t\t   PLoS Genetics, 4(10):e1000212, 2008." << endl;
	cout << "\t\t   doi:10.1371/journal.pgen.1000212" << endl;
        return 0;
    }

    ifstream raw_mapping_data_file(argv[1]);
    string tmp_str;
    string population_type;

    raw_mapping_data_file >> tmp_str;
    if (tmp_str != "population_type")
    {
        cout << "ERROR, the input file is invalid" << endl;
        assert(false); // crash the program on error
        return -1;
    }
    raw_mapping_data_file >> population_type;
    raw_mapping_data_file.close();

    genetic_map* barley;
    if (population_type == "DH") {
        barley = new genetic_map_DH();
    } else {
        barley = new genetic_map_RIL();    
    }
    
    
    barley->read_raw_mapping_data(argv[1]);

    /*The algorithm parameter is provided*/
    barley->generate_map();
    
    ofstream output_file(argv[2]);
    barley->write_output(output_file);
    output_file.close();
    
    delete barley;
    return 0;
}

void print_vector(vector<int> tmp)
{
    for (int ii = 0 ; ii < tmp.size(); ii++)
    {
        cout << ii << ',' ;
    }
    cout << endl;
}
