#include <map>
using namespace std;


typedef enum {SYS_NOMINAL,
	SYS_PU_UP, SYS_PU_DOWN,
	SYS_TOP_PT} systematic_shift;

systematic_shift allSystematics[] = {SYS_NOMINAL,
	SYS_PU_UP, SYS_PU_DOWN,
	SYS_TOP_PT};

// C++11 feature:
map<systematic_shift, const char*> systematic_shift_names = {{SYS_NOMINAL, "NOMINAL"},
	{SYS_PU_UP, "PU_UP"}, {SYS_PU_DOWN, "PU_DOWN"},
	{SYS_TOP_PT, "TOP_PT"}};

/*
const char * systematic_shift_name[] = {
    "NOMINAL",
    "PU_UP",
    "PU_DOWN",
};
*/

