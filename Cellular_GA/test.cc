
#include <sstream>
#include <iostream>
#include <fstream>

#include <cstdlib>
#include <ctime>

#include <map>

using namespace std;



int main(void)
{
    map<unsigned char, int> test;

    test['a']++;

    cout << test['a'];
}
