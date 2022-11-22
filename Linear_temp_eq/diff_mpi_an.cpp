#include<iostream>
#include<fstream>
#include<string>
#include<stdlib.h>

using namespace std;
bool tokenize(string s1, string s2, string del = " ")
{
    int start, end = -1*del.size();
    bool check = 1;
    do{
        start = end + del.size();
        end = s1.find(del, start);
        if (atof((s1.substr(start, end - start)).c_str()) and atof((s2.substr(start, end - start)).c_str()))
            if (atof((s1.substr(start, end - start)).c_str()) - atof((s2.substr(start, end - start)).c_str()) > 0){ 
                cout << "Difference is found!" << endl;
                check = 0;
                break;
            }
    }while (end != -1);
    return check;
}

int main() {
    ifstream analit("file_compare_with_mpi.dat");
    ifstream num("file_mpi.dat");
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