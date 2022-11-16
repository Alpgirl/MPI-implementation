#include<iostream>
#include<fstream>
#include<string>
#include<stdlib.h>

using namespace std;
void tokenize(string s1, string s2, string del = " ")
{
    int start, end = -1*del.size();
    do{
        start = end + del.size();
        end = s1.find(del, start);
        if (atof((s1.substr(start, end - start)).c_str()) and atof((s2.substr(start, end - start)).c_str()))
            cout << atof((s1.substr(start, end - start)).c_str()) - atof((s2.substr(start, end - start)).c_str()) << " ";
    }while (end != -1);
}

int main() {
    ifstream analit("file_compare_with_mpi.dat");
    ifstream num("file_mpi.dat");
    string line_an, line_num, del = " ";
    double diff;
    while (getline(analit,line_an) and getline(num, line_num)){
       // while ((pos = line_an.find(del)) != '\n') {
            tokenize(line_an, line_num, " ");
            cout << endl;
            //diff = atof((line_an.substr(0, pos)).c_str())-atof((line_num.substr(0, pos)).c_str());
            //cout << atof((line_an.substr(0, pos)).c_str()) << " ";
            //line_an.erase(0, pos + del.length());
           // line_num.erase(0, pos + del.length());
        //}
    }
    analit.close();
    num.close();
    return 0;
}