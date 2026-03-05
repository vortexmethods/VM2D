/*--------------------------------*- VM2D -*-----------------*---------------*\
| ##  ## ##   ##  ####  #####   |                            | Version 1.14   |
| ##  ## ### ### ##  ## ##  ##  |  VM2D: Vortex Method       | 2026/03/06     |
| ##  ## ## # ##    ##  ##  ##  |  for 2D Flow Simulation    *----------------*
|  ####  ##   ##   ##   ##  ##  |  Open Source Code                           |
|   ##   ##   ## ###### #####   |  https://www.github.com/vortexmethods/VM2D  |
|                                                                             |
| Copyright (C) 2017-2026 I. Marchevsky, K. Sokol, E. Ryatina, A. Kolganova   |
*-----------------------------------------------------------------------------*
| File name: TimesGen.h                                                       |
| Info: Source code of VM2D                                                   |
|                                                                             |
| This file is part of VM2D.                                                  |
| VM2D is free software: you can redistribute it and/or modify it             |
| under the terms of the GNU General Public License as published by           |
| the Free Software Foundation, either version 3 of the License, or           |
| (at your option) any later version.                                         |
|                                                                             |
| VM2D is distributed in the hope that it will be useful, but WITHOUT         |
| ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or       |
| FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License       |
| for more details.                                                           |
|                                                                             |
| You should have received a copy of the GNU General Public License           |
| along with VM2D.  If not, see <http://www.gnu.org/licenses/>.               |
\*---------------------------------------------------------------------------*/


/*!
\file
\brief Заголовочный файл с описанием класса TimesGen
\author Марчевский Илья Константинович
\author Сокол Ксения Сергеевна
\author Рятина Евгения Павловна
\author Колганова Александра Олеговна
\Version 1.14
\date 6 марта 2026 г.
*/ 

#ifndef TIMESGEN_H
#define TIMESGEN_H

#include <chrono>
#include <map>
#include <memory>
#include "defs.h"

namespace VMlib
{   
    class WorldGen;
    
    /*!
    \brief Класс засекания времени
    \author Марчевский Илья Константинович
    \Version 1.14
    \date 6 марта 2026 г.
    */
    class vmTimer
    {
        /// Имя счетчика времени
        const std::string label;

        using Clock = std::chrono::high_resolution_clock;

        /// Признак того, что счетчик работает
        mutable bool active = false;

        /// Хранит накопденную от предыдущих запусков продолжительность работы счетчика
        mutable std::chrono::duration<float> duration_;

        /// Отметка последнего включения счетчика
        mutable Clock::time_point start_ = Clock::now();

        /// Отметка последнего выключения счетчика
        mutable Clock::time_point stop_ = Clock::now();

        // Запрет на копирование счетчика
        vmTimer(const vmTimer&) = delete;
        vmTimer& operator=(const vmTimer&) = delete;

    public:
        //using ns = std::chrono::nanoseconds;
        //using mks = std::chrono::microseconds;
        //using ms = std::chrono::milliseconds;
        //using s = std::chrono::seconds;

        /// Основная единица измерения времени --- миллисекунды
        using ms = std::chrono::duration<float, std::chrono::milliseconds::period>;
        
        /// Дополнительная единица измерения времени --- секунды
        using s = std::chrono::duration<float, std::chrono::seconds::period>;

        /// \brief Конструктор, принимающий на вход имя счетчика
        ///
        /// Создает счетчик, засекающий время в миллисекундах
        /// \param[in] timerLabel константная ссылка на строку --- имя счетчика
        vmTimer(const std::string& timerLabel = "") : label(timerLabel) { reset(); }
        ~vmTimer() = default;
        
        /// Сброс счетчика времени
        const vmTimer& reset() const
        {
            duration_ = std::chrono::seconds(0);
            active = false;
            return *this;
        }

        /// Запуск (первый или повторный) счетчика времени
        const vmTimer& start() const
        {
            if (!active)
            {
                start_ = Clock::now();
                active = true;
            }
            else
            {
                std::cout << "Timer " << label << " was not stopped before being start!" << std::endl;
                exit(-15);
            }
            return *this;
        }

        /// Останов работающего счетчика времени
        const vmTimer& stop() const
        {
            if (active)
            {
                stop_ = Clock::now();
                duration_ += stop_ - start_;
                active = false;
            }
            else
            {
                std::cout << "Timer " << label << " was not started before being stopped!" << std::endl;
                exit(-15);
            }
            return *this;
        }


        template<typename T = ms>
        double duration() const
        {
            return std::chrono::duration_cast<T>(duration_).count();
        }
    };


    /*!
    \brief Класс для сбора статистики времени исполнения основных шагов алгоритма и вывода ее в файл
    \author Марчевский Илья Константинович
    \Version 1.14
    \date 6 марта 2026 г.
    */
    class TimersGen    
    { 
    private:
        /// Константная ссылка на решаемую задачу
        const WorldGen& W;

    protected:
        /// Список имен счетчиков
        std::vector<std::string> timerLabelList;

        /// Ассоциативный массив { имя, счетчик }
        std::map<std::string, std::unique_ptr<vmTimer>> timer;    
    
    public:
        /// Конструктор
        TimersGen(const WorldGen& W_, std::vector<std::string> labels);

        /// Запуск счетчика
        void start(const std::string& timerLabel);

        /// Останов счетчика
        void stop(const std::string& timerLabel);

        /// Сброс всех счетчиков
        void resetAll();

        /// Вывод счетчика всего шага в секундах
        template<typename T = vmTimer::s>
        double durationStep() const
        {
            return (timer.at("Step"))->duration<T>();
        }

        /// Формирование заголовка файла временной статистики
        void GenerateStatHeader();

        /// Формирование очередной строки файла временной статистики
        void GenerateStatString(size_t stepNo, double curTime, size_t N);
    };




	class TimesGen
	{
	protected:
		/// Обнуление одного временного периода 
		/// \param[out] period промежуток времени, начало и конец которого будут обнулены
		static void ToZeroPeriod(timePeriod& period)
		{
			period.first = 0;
			period.second = 0;
		}//ToZero(...)

	public:
		/// Конструктор
		TimesGen() {};

		/// Деструктор
		virtual ~TimesGen() {};


		/// Генерация заголовка файла временной статистики
		virtual void GenerateStatHeader() const = 0;

		/// Сохранение строки со статистикой в файл временной статистики
		virtual void GenerateStatString() const = 0;
		
		/// Обнуление состояния временной статистики
		virtual void ToZero() = 0;

		/// Вычисление разницы во времени для пары засечек в секундах
		/// \param[in] t константная ссылка на пару засечек времени
		/// \return разницу в секундах
		static double dT(const timePeriod& t)
		{
			return (t.second - t.first);
		}//dT(...)
	};

}//namespace VMlib

#endif