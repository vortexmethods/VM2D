#include "CommentsParser.h"

char CommentsParser::processSymbol(char ch)
{
	switch (ch)
	{
	case '\n':
		return ' ';

	case ';':
		return '\n';

	default:
		return ch;
	}
}

void CommentsParser::normalParser(char ch, std::string& str)
{
	switch (ch)
	{
	case '\"':
		currentParser = &CommentsParser::inStringParser;
		break;

	case '\'':
		currentParser = &CommentsParser::inCharParser;
		break;

	case '/':
		currentParser = &CommentsParser::afterSlashParser;
		return;
	}
	str.push_back(processSymbol(ch));
}


void CommentsParser::inStringParser(char ch, std::string& str)
{
	str.push_back(ch);
	if (ch == '\"' && !inStringLastEscape)
		currentParser = &CommentsParser::normalParser;
	inStringLastEscape = (ch == '\\') && !inStringLastEscape;
}


void CommentsParser::inCharParser(char ch, std::string& str)
{
	str.push_back(ch);
	if (ch == '\'' && !inCharLastEscape)
		currentParser = &CommentsParser::normalParser;
	inCharLastEscape = (ch == '\\') && !inCharLastEscape;
}

void CommentsParser::afterSlashParser(char ch, std::string& str)
{
	switch (ch)
	{
	case '/':
		currentParser = &CommentsParser::inInlineCommentParser;
		return;

	case '*':
		currentParser = &CommentsParser::inMultilineCommentParser;
		return;

	default:
		str.push_back(processSymbol(ch));
	}
}


void CommentsParser::inInlineCommentParser(char ch, std::string& str)
{
	if (ch == '\n' && !inInlineLastEscape)
	{
		str.push_back(processSymbol(ch));
		currentParser = &CommentsParser::normalParser;
	}
	inInlineLastEscape = (ch == '\\') && !inInlineLastEscape;
}


void CommentsParser::inMultilineCommentParser(char ch, std::string& str)
{
	if (ch == '/' && inMultilineLastStar)
	{
		str.push_back(' ');
		currentParser = &CommentsParser::normalParser;
	}
	inMultilineLastStar = (ch == '*');
}