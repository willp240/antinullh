#ifndef __ANTINUFIT__SystFactory__
#define __ANTINUFIT__SystFactory__
#include <string>

class Systematic;
namespace antinufit
{
	class SystFactory
	{
	public:
		static Systematic *New(const std::string &name,
													 const std::string &type_,
													 double value);
	};
}
#endif
