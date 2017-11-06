/*!
\file
\brief Заголовочный файл с описанием класса StreamParser
\author Марчевский Илья Константинович
\author Кузьмина Ксения Сергеевна
\author Рятина Евгения Павловна
\version 1.0
\date 1 декабря 2017 г.
*/

#ifndef STREAMPARSER_H
#define STREAMPARSER_H

//#define VARNAME(var) #var

#include <unordered_map>
#include <sstream>
#include <algorithm>

#include "defs.h"
#include "Vortex2D.h"

/*!
\brief Класс, позволяющий выполнять разбор файлов и строк с настройками и параметрами
\author Марчевский Илья Константинович
\author Кузьмина Ксения Сергеевна
\author Рятина Евгения Павловна
\version 1.0
\date 1 декабря 2017 г.
*/
class StreamParser
{
private:
	/// Указатель на поток вывода логов
	static std::ostream* pinfo;

	/// Указатель на поток вывода ошибок
	static std::ostream* perr;

	/// \brief Данные считываемые из основного файла-паспорта расчета
	///
	/// Неупорядоченный ассоциативный контейнер, хранящий данные из паспорта расчета
	/// \n По умолчанию данные берутся из файла passport.txt, но это имя может быть изменено ключом passport pspFile в строке параметров задачи 
	std::unordered_map<std::string, std::vector<std::string>> database;
	
	/// \brief Данные считываемые из списка умолчаний
	///
	/// Неупорядоченный ассоциативный контейнер, хранящий данные из файла со списком значений по умолчанию
	/// \n По умолчанию данные берутся из файла defaults.txt
	std::unordered_map<std::string, std::vector<std::string>> defaults;

	/// \brief Данные считываемые из параметров конкретной задачи
	///
	/// Неупорядоченный ассоциативный контейнер, хранящий данные из строки --- перечня параметров конкретной задачи в файле tasks.txt. 
	/// \n Здесь же определяются переменные, которые затем могут использоваться методом разыменовывания в файлах-паспортах конкретных задач
	std::unordered_map<std::string, std::vector<std::string>> vars;

	/// \brief Данные считываемые из перечня параметров-переключателей
	///
	/// Неупорядоченный ассоциативный контейнер, хранящий значения параметров-переключателей 
	/// \n По умолчанию данные берутся из файла switchers.txt. 	
	std::unordered_map<std::string, std::vector<std::string>> switchers;

	/// \brief Парсинг заданного потока 
	///
	/// Разбирает поток вида key = value   или   key = { value1, value2, ..., valueN }
	/// \n Возвращает базу данных и список ключей в порядке появления, указывается тип скобок, по умолчанию --- круглые
	/// \param[in] stream ссылка на поток, который требуется распарсить
	/// \param[out] database ссылка на заполняемую базу данных (unordered_map)
	/// \param[out] vekKey ссылка на список ключей в порядке появления
	/// \param[in] openBracket тип открывающейся скобки (по умолчанию --- "(" )
	/// \param[in] closeBracket тип закрывающейся скобки (по умолчанию --- ")" )
	void ParseStream(std::istream& stream, std::unordered_map<std::string, std::vector<std::string>>& database, std::vector<std::string>& vekKey, char openBracket = '(', char closeBracket = ')');

public:
	/// \brief Конструктор, принимающий четыре потока
	/// 
	/// Заполняет 4 базы данных --- database, defaults, switchers, vars
	/// \n Используется для установки основных параметров расчета
	///
	/// \param[in] mainStream ссылка на основной поток параметров
	/// \param[in] defaultsStream ссылка на поток по умолчанию
	/// \param[in] switchersStream ссылка на поток параметров-переключателей
	/// \param[in] varsStream ссылка на поток переменных
	StreamParser(std::istream& mainStream, std::istream& defaultsStream, std::istream& switchersStream, std::istream& varsStream);
		
	/// \brief Конструктор, принимающий два потока
	/// 
	/// Заполняет 2 базы данных --- database, defaults; базы данных switchers и vars остаются пустыми
	/// \n Используется для первоначального считывания параметров задач из списка задач и формирования очереди
	///
	/// \param[in] mainStream ссылка на основной поток параметров
	/// \param[in] defaultsStream ссылка на поток по умолчанию
	StreamParser(std::istream& mainStream, std::istream& defaultsStream);

	/// \brief Конструктор, принимающий два потока и возвращающий список ключей
	/// 
	/// Заполняет 2 базы данных --- database, defaults; базы данных switchers и vars остаются пустыми; 
	/// возвращает список ключей в порядке их появления
	/// \n Используется для обработки списка задач
	///
	/// \param[in] mainStream ссылка на основной поток параметров
	/// \param[in] defaultsStream ссылка на поток по умолчанию
	/// \param[out] vekKeys ссылка на список ключей в порядке их появления
	StreamParser(std::istream& mainStream, std::istream& defaultsStream, std::vector<std::string>& vekKeys);
		
	/// \brief Конструктор, принимающий один поток
	/// 
	/// Заполняет 1 базу данных --- database; базы данных defaults, switchers и vars остаются пустыми
	/// \n Используется для считывания профиля и следа
	///
	/// \param[in] mainStream ссылка на основной поток параметров
	StreamParser(std::istream& mainStream);

	/// Деструктор
	~StreamParser() { };
	
	/// \brief Pазбор строки, содержащей запятые, на отдельные строки
	///
	/// Запятые внутри парных скобок игнорируются и не воспринимаются как разделители
	/// \param[in] line разбираемая строка
	/// \param[in] openBracket тип открывающейся скобки (по умолчанию --- "(" )
	/// \param[in] closeBracket тип закрывающейся скобки (по умолчанию --- ")" )
	/// \return список (вектор) из отдельных строк
	static std::vector<std::string> StringToVector(std::string line, char openBracket = '(', char closeBracket = ')');
	
	/// \brief Объединение вектора (списка) из строк в одну строку
	///
	/// После каждой строки вставляется символ конца строки
	/// \param[in] _vecString константная ссылка на вектор из строк
	/// \return сплошную строку с символами конца строки между подстроками
	static std::string VectorStringToString(const std::vector<std::string>& _vecString);
	
	/// \brief Разбор строки на пару ключ-значение
	///
	/// Разбирает строку вида xxx(yyy) на пару подстрок xxx и yyy
	/// \param[in] line разбираемая строка вида xxx(ууу)
	/// \return пара строк <xxx, yyy>
	static std::pair<std::string, std::string> SplitString(std::string line);

	/// \brief Установка значения параметра по умолчанию
    ///
	/// Устанавливает параметру значение по умолчанию, если оно задано в виде константной ссылки
	/// \tparam T тип параметра, значение которого устанавливается
	/// \param[in] name константная ссылка на строку --- имя устанавливаемого параметра (необходимо только для отображения в лог)
	/// \param[in] res ссылка на задаваемый параметр
	/// \param[in] defValue указатель на константу --- значение по умолчанию
	template <typename T>
	void SetDefault(const std::string& name, T& res, const T* defValue) const
	{
		if (defValue != nullptr)
		{
			res = *defValue;
			*pinfo << "parser info: parameter <" << name << " = " << res << "> set as default" << std::endl;
		}
		else
		{
			*perr << "parser error: PARAMETER " << name << " NOT FOUND " << std::endl;
			exit(-1);
		}
	}

	/// \brief Считывание вектора из двумерных точек из базы данных
	/// 
	/// Переопределение метода get() для ситывания вектора из точек (Point2D) из базы данных
	/// \param[in] name константная ссылка на строку --- имя считываемого параметра
	/// \param[out] res ссылка на значение, считываемое из базы данных
	/// \param[in] defValue указатель на константу --- значение по умолчанию (по умолчанию nullptr)
	void get(const std::string& name, std::vector<Point2D>& res, const std::vector<Point2D>* defValue = nullptr) const
	{
		res.resize(0);
		std::unordered_map<std::string, std::vector<std::string>>::const_iterator search_it;

		//поиск ключа в базе и, если нет, то в базе умолчаний
		if (((search_it = database.find(name)) != database.end()) ||
			((search_it = defaults.find(name)) != defaults.end()))
		{
			for (size_t i = 0; i < (search_it->second).size(); ++i)
			{
				std::string s = (search_it->second)[i];
				//проверка на подстановку
				if (s[0] != '$')
				{
					//если не подстановка
					//stringstream ss(s);
					Point2D elem;
					//ss >> elem;
					int pos = s.find(',', 1);
					std::stringstream(s.substr(1, pos)) >> elem[0];
					std::stringstream(s.substr(pos+1, s.length()-pos)) >> elem[1];
					res.push_back(elem);
				}
				else
				{
					//если подстановка
					*perr << "parser error: VECTOR SUBSTITUTIONS ARE NOT PERMITED (FOR TPOINT2D TOO)!" << std::endl;
					exit(-1);
				}
			}
		}
		//else
		//	SetDefault(name, res, defValue);
	};
		
	/// \brief Считывание вектора из вихрей из базы данных
	/// 
	/// Переопределение метода get() для ситывания вектора из вихрей (Vortex2D) из базы данных
	/// \param[in] name константная ссылка на строку --- имя считываемого параметра
	/// \param[out] res ссылка на вектор из вихрей, считываемый из базы данных
	/// \param[in] defValue указатель на константу --- значение по умолчанию (по умолчанию nullptr)
	void get(const std::string& name, std::vector<Vortex2D>& res, const std::vector<Vortex2D>* defValue = nullptr) const
	{
		res.resize(0);
		std::unordered_map<std::string, std::vector<std::string>>::const_iterator search_it;

		//поиск ключа в базе и, если нет, то в базе умолчаний
		if (((search_it = database.find(name)) != database.end()) ||
			((search_it = defaults.find(name)) != defaults.end()))
		{
			for (size_t i = 0; i < (search_it->second).size(); ++i)
			{
				std::string s = (search_it->second)[i];
				//проверка на подстановку
				if (s[0] != '$')
				{
					//если не подстановка
					Point2D r;
					double g;
					int pos1 = s.find(',', 1);
					int pos2 = s.find(',', pos1 + 1);
					std::stringstream(s.substr(1, pos1)) >> r[0];
					std::stringstream(s.substr(pos1 + 1, s.length() - pos1)) >> r[1];
					std::stringstream(s.substr(pos2 + 1, s.length() - pos2)) >> g;

					Vortex2D elem(r, g);

					res.push_back(elem);
				}
				else
				{
					//если подстановка
					*perr << "parser error: VECTOR SUBSTITUTIONS ARE NOT PERMITED (FOR TVORTEX2D TOO)!" << std::endl;
					exit(-1);
				}
			}
		}
		//else
		//	SetDefault(name, res, defValue);
	};
	
	/// \brief Считывание вектора из простых типов из базы данных
	/// 
	/// Шаблонный метод get() для считывания вектора из простых типов из базы данных
	/// \tparam T тип считываемых данных
	/// \param[in] name константная ссылка на строку --- имя считываемого параметра
	/// \param[out] res ссылка на вектор из данных, считываемый из базы данных
	/// \param[in] defValue указатель на константу --- значение по умолчанию (по умолчанию nullptr)
	template <typename T>
	void get(const std::string& name, std::vector<T>& res, const std::vector<T>* defValue = nullptr) const
	{
		res.resize(0);
		std::unordered_map<std::string, std::vector<std::string>>::const_iterator search_it;

		//поиск ключа в базе и, если нет, то в базе умолчаний
		if ( ((search_it = database.find(name)) != database.end()) ||
			 ((search_it = defaults.find(name)) != defaults.end()))
		{
			for (size_t i = 0; i < (search_it->second).size(); ++i)
			{
				std::string s = (search_it->second)[i];
				//проверка на подстановку
				if (s[0] != '$')
				{
					//если не подстановка
					std::stringstream ss(s);
					T elem;
					ss >> elem;
					res.push_back(elem);
				}
				else
				{
					//если подстановка
					*perr << "parser error: VECTOR SUBSTITUTIONS ARE NOT PERMITED!" << std::endl;
					exit(-1);
				}
			}
		}
		else
			SetDefault(name, res, defValue);
	};	
	
	/// \brief Считывание скаляра из базы данных
	/// 
	/// Шаблонный метод get() для считывания скаляра из базы данных
	/// \tparam T тип считываемых данных
	/// \param[in] name константная ссылка на строку --- имя считываемого параметра
	/// \param[out] res ссылка данные, считываемые из базы данных
	/// \param[in] defValue указатель на константу --- значение по умолчанию (по умолчанию nullptr)
	template <typename T>
	void get(const std::string& name, T& res, const T* defValue = nullptr) const
	{
		std::unordered_map<std::string, std::vector<std::string>>::const_iterator search_it;

		//поиск ключа в базе и, если нет, то в базе умолчаний
		if ( (( (search_it = database.find(name))  != database.end()) && ((search_it->second).size() > 0)) ||
			 (( (search_it = defaults.find(name))  != defaults.end()) && ((search_it->second).size() > 0))  )
		{
			if ((search_it->second).size() == 1)
			{
				std::string s = (search_it->second)[0];
				
				//проверка на значение-переключатель
				if ( ((search_it = switchers.find(s)) != switchers.end()) && ((search_it->second).size() > 0) )
				{
					//если переключатель
					std::stringstream ss(search_it->second[0]);
					T elem;
					ss >> elem;
					res = elem;
				}
				
				//проверка на подстановку
				else if (s[0] != '$')
				{
					//если не подстановка
					std::stringstream ss(s);
					T elem;
					ss >> elem;
					res = elem;
				}
				else
				{
					//если подстановка
					std::string svar = s.substr(1, s.length() - 1);
					auto search_var = vars.find(svar);
					if ((search_var != vars.end()) && ((search_var->second).size() > 0))
					{
						std::stringstream ss(search_var->second[0]);
						T elem;
						ss >> elem;
						res = elem;
					}
					else
					{
						*perr << "parser error: SUBSTITUTION $" << svar << " IS UNDEFINED" << std::endl;
						exit(1);
					}
				}
			}
			else
			{
				*perr << "parser error: PARAMETER " << name << " IS VECTOR " << std::endl;
				exit(-1);
			}
		}
		else
			SetDefault(name, res, defValue);
	};

	/// \brief Считывание точки из базы данных
	/// 
	/// Переопределение метода get() для ситывания точки (Point2D) из базы данных
	/// \param[in] name константная ссылка на строку --- имя считываемого параметра
	/// \param[out] res ссылка на точку, считываемую из базы данных
	/// \param[in] defValue указатель на константу --- значение по умолчанию (по умолчанию nullptr)
	void get(const std::string& name, Point2D& res, const Point2D* defValue = nullptr) const
	{
		std::unordered_map<std::string, std::vector<std::string>>::const_iterator search_it;

		//поиск ключа в базе и, если нет, то в базе умолчаний
		if ( ((search_it = database.find(name)) != database.end()) ||
			 ((search_it = defaults.find(name)) != defaults.end()) )
		{
			if ((search_it->second).size() == res.size())
			{
				for (size_t i = 0; i < (search_it->second).size(); ++i)
				{
					std::string s = (search_it->second)[i];
					//проверка на подстановку
					if (s[0] != '$')
					{
						//если не подстановка
						std::stringstream ss(s);
						double elem;
						ss >> elem;
						res[i] = elem;
					}
					else
					{
						//если подстановка
						std::string svar = s.substr(1, s.length() - 1);
						auto search_var = vars.find(svar);
						if ((search_var != vars.end()) && ((search_var->second).size() > 0))
						{
							std::stringstream ss(search_var->second[0]);
							double elem;
							ss >> elem;
							res[i] = elem;
						}
						else
						{
							*perr << "parser error: SUBSTITUTION $" << svar << " IS UNDEFINED" << std::endl;
							exit(-1);
						}
					}
				}
			}
			else
			{
				*perr << "parser error: PARAMETER " << name << " LENGTH DIFFERS FROM 2" << std::endl;
				exit(-1);
			}

		}
		else
			SetDefault(name, res, defValue);
	}


};

#endif