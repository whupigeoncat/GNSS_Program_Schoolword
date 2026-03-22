#include <iostream>
#include <iomanip>
#include <cmath>
#include "RTK_Structs.h"

using namespace std;

static bool NearlyEqual(double a, double b, double tol)
{
	return fabs(a - b) <= tol;
}

int main()
{
    std::string filepath = R"(C:\Users\asus\Documents\Program\Program_GNSS\RTK_CSH\Data\oem719-202603111200.bin)";
    if (!DataReceive(filepath))
    {
        cout << "false";
        return 1;
    }
    return 0;
}