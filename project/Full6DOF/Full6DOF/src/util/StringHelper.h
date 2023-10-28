#pragma once

#include <vector>
#include <string>

namespace strings
{
	std::vector<std::string> split(const std::string& input, char deliminator = ' ')
	{
		std::vector<std::string> tokens;
		int token_start = 0;
		int token_end = 0;
		while (token_start < input.size())
		{
			while (input[token_start] == deliminator)
			{
				token_start++;
			}
			token_end = token_start;
			while (input[token_end] != deliminator)
			{
				token_end++;
			}
			tokens.push_back(input.substr(token_start, token_end - token_start));
			token_start = token_end;
		}
		return tokens;
	}

}