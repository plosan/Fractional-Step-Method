#include <iostream>
#include <fstream>
#include <string>
#include <sstream>

using namespace std;


int main() {

    string name1 = "benchmark_Re100_129_129/p_129_129_100.txt";
    string name2 = "test/p_129_129_100.txt";

    cout << "File 1: " << name1 << endl;
    cout << "File 2: " << name2 << endl;

    ifstream f1, f2;
    f1.open(name1);
    if(!f1.is_open()) {
        cout << "Error opening file 1" << endl;
        f1.close();
    }

    f2.open(name2);
    if(!f2.is_open()) {
        cout << "Error opening file 2" << endl;
        f2.close();
    }

    string line1;
    bool found = false;
    int count = 0;
    while(getline(f1, line1) && !found) {

        count++;
        string line2;
        getline(f2, line2);

        int comp = line1.compare(line2);
        // cout << comp << endl;

        if(comp != 0) {
            found = true;
            cout << count << " File1 " << line1 << endl;
            cout << count << " File2 " << line2 << endl;
        }

    }

    f1.close();
    f2.close();



}
