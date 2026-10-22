#pragma once

#include <string>
#include <vector>

namespace antinufit
{
  std::vector<std::string> SplitString(const std::string &input, char delimiter);
  std::string stripQuoteMarks(std::string value);
}
