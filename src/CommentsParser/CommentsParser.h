/*!
\file
\brief Заголовочный файл с описанием класса CommentsParser
\author Марчевский Илья Константинович
\author Кузьмина Ксения Сергеевна
\author Рятина Евгения Павловна
\version 1.0
\date 1 декабря 2017 г.
*/

#ifndef COMMENTSPARSER_H
#define COMMENTSPARSER_H

#include <string>
#include <sstream>
#include <algorithm>

class CommentsParser
{
private:
	void(CommentsParser::*currentParser)(char ch, std::string& str);
	
	void normalParser(char ch, std::string& str);
	void inStringParser(char ch, std::string& str);
	void inCharParser(char ch, std::string& str);
	void afterSlashParser(char ch, std::string& str);
	void inInlineCommentParser(char ch, std::string& str);
	void inMultilineCommentParser(char ch, std::string& str);

	bool inStringLastEscape;
	bool inCharLastEscape;
	bool inInlineLastEscape;
	bool inMultilineLastStar;

	char CommentsParser::processSymbol(char ch);

public:
	std::string intermediateOutput;
	std::stringstream resultStream;

	CommentsParser(const std::string& fileName)
		:
		inStringLastEscape(false), 
		inCharLastEscape(false), 
		inInlineLastEscape(false),
		inMultilineLastStar(false),
		currentParser(&CommentsParser::normalParser)
	{
		FILE * inputFile;
		fopen_s(&inputFile, fileName.c_str(), "r");

		int symbol;
		
		while ((symbol = fgetc(inputFile)) != EOF)
			(this->*currentParser)(symbol, intermediateOutput);

		fclose(inputFile);

		std::stringstream ss(intermediateOutput);
		
		std::string readline;
		
		while (ss.good())
		{
			getline(ss, readline);
		
			readline.erase(remove(readline.begin(), readline.end(), ' '), readline.end());
			readline.erase(remove(readline.begin(), readline.end(), '\t'), readline.end());
	
			resultStream << readline << std::endl;
		} 
	};

	CommentsParser(){ };

	~CommentsParser(){ };
};

#endif