#include <StringUtils.hh>

#include <sstream>

namespace antinufit
{
  std::vector<std::string> SplitString(const std::string &input, char delimiter)
  {
    std::vector<std::string> parts;
    std::stringstream stream(input);
    std::string part;
    while (std::getline(stream, part, delimiter))
      parts.push_back(part);
    return parts;
  }

  std::string stripQuoteMarks(std::string value)
  {
    if (value.empty())
      return value;
    size_t first = 0;
    size_t last = value.size() - 1;
    if (value[first] == '"' || value[first] == '\'')
      ++first;
    if (last >= first && (value[last] == '"' || value[last] == '\''))
      --last;
    return first > last ? std::string() : value.substr(first, last - first + 1);
  }
}
