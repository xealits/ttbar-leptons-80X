#include <map>
using namespace std;


typedef enum {SYS_NOMINAL,
	SYS_PU_UP, SYS_PU_DOWN,
	SYS_TOP_PT,
	SYS_JER_UP,
	SYS_JER_DOWN} systematic_shift;

systematic_shift weightSystematics[] = {SYS_NOMINAL,
	SYS_PU_UP, SYS_PU_DOWN,
	SYS_TOP_PT};

systematic_shift jetSystematics[] = {SYS_NOMINAL,
	SYS_JER_UP,
	SYS_JER_DOWN};

systematic_shift allSystematics[] = {SYS_NOMINAL,
	SYS_PU_UP, SYS_PU_DOWN,
	SYS_TOP_PT,
	SYS_JER_UP,
	SYS_JER_DOWN};

// C++11 feature:
map<systematic_shift, const char*> systematic_shift_names = {{SYS_NOMINAL, "NOMINAL"},
	{SYS_PU_UP, "PU_UP"}, {SYS_PU_DOWN, "PU_DOWN"},
	{SYS_TOP_PT, "TOP_PT"},
	{SYS_JER_UP,   "JER_UP"},
	{SYS_JER_DOWN, "JER_DOWN"} };

/*
const char * systematic_shift_name[] = {
    "NOMINAL",
    "PU_UP",
    "PU_DOWN",
};
*/

