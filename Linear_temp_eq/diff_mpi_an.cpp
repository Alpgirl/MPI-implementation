#include<iostream>
#include<fstream>
#include<string>
#include<stdlib.h>
#include<math.h>

using namespace std;
bool tokenize(string s1, string s2, string del = " ")
{
    int start1, start2, end1 = -1*del.size(), end2 = -1*del.size();
    bool check = 1;
    double tmp1, tmp2;
    do{
        start1 = end1 + del.size();
        start2 = end2 + del.size();
        end1 = s1.find(del, start1);
        end2 = s2.find(del, start2);
        tmp1 = atof((s1.substr(start1, end1 - start1)).c_str());
        tmp2 = atof((s2.substr(start2, end2 - start2)).c_str());
        if (tmp1 and tmp2)
            if (abs(tmp1 - tmp2) > 0){ 
                cout << "Difference is found!" << endl;
                cout << s1.substr(start1, end1 - start1) << " " << s2.substr(start2, end2 - start2) << " " << start1 << " " << start2 << endl;
                check = 0;
                break;
            }
    }while (end1 != -1 or end2 != -1);
    return check;
}

int main() {
    //ifstream analit("file_compare_with_mpi.dat");
    ifstream analit("file_analytical_compare_mpi.dat");
    //ifstream num("file_mpi.dat");
    ifstream num("file_mpi_analyt.dat");
    string line_an, line_num, del = " ";
    double diff;
    bool check = 1;
    int cnt = 0;
    while (getline(analit,line_an) and getline(num, line_num)){
            check = tokenize(line_an, line_num, " ");
            if (check == 0) cnt++;
    }
    if (cnt == 0) cout << "Good job!" << endl;
    analit.close();
    num.close();
    return 0;
}