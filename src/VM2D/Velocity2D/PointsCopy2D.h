/*--------------------------------*- VM2D -*-----------------*---------------*\
| ##  ## ##   ##  ####  #####   |                            | Version 1.9    |
| ##  ## ### ### ##  ## ##  ##  |  VM2D: Vortex Method       | 2020/07/22     |
| ##  ## ## # ##    ##  ##  ##  |  for 2D Flow Simulation    *----------------*
|  ####  ##   ##   ##   ##  ##  |  Open Source Code                           |
|   ##   ##   ## ###### #####   |  https://www.github.com/vortexmethods/VM2D  |
|                                                                             |
| Copyright (C) 2017-2020 Ilia Marchevsky, Kseniia Kuzmina, Evgeniya Ryatina  |
*-----------------------------------------------------------------------------*
| File name: PointsCopy2D.h                                                   |
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
\brief Заголовочный файл с описанием класса PointsCopy
\author Марчевский Илья Константинович
\author Кузьмина Ксения Сергеевна
\author Рятина Евгения Павловна
\version 1.9
\date 22 июля 2020 г.
*/


#ifndef POINTSCOPY_H_
#define POINTSCOPY_H_

#include "Point2D.h"
#include "Vortex2D.h"
#include "defs.h"

namespace PointsCopy_iterator
{
	struct all_reference;
	struct all_copy;
	struct all_iterator;
}
namespace VM2D
{
	/*!
	\brief Класс, опеделяющий копию для сортировки структуры данных следа (пелены), источников или "виртуальных" вихрей или источников
	\author Марчевский Илья Константинович
	\author Кузьмина Ксения Сергеевна
	\author Рятина Евгения Павловна
	\version 1.9
	\date 22 июля 2020 г.
	*/
class PointsCopy
{
	public:
		/// Список из самих объектов, приведенных к типу Vortex2D 
		std::vector<Vortex2D> vtx;
		
		///\brief Список из пар значений - номер рассматриваемого профиля и номер панели на нем
		/// Для wake, sourceWake или wakeVP не задается 
		std::vector<std::pair<int, int>> aflPnl;

		/// Список из номеров элементов из vtx в исходном массиве точек 
		std::vector<int> num;

		PointType type;

		/// Ниже - несортируемые списки:
		std::vector<Point2D> velo;
		std::vector<double> domainRadius;

		void push_back(Vortex2D _vtx, std::pair<int, int> _aflPnl, int _num);
		void reserve(size_t size);
		void emplace_back(Vortex2D _vtx, std::pair<int, int> _aflPnl, int _num, PointType _type, Point2D _velo = { 0.0, 0.0 }, double _domainRadius = 0.0);
		void clear();

		using iterator = PointsCopy_iterator::all_iterator;
		iterator begin();
		iterator end();

		size_t size() const;

	};
}

namespace PointsCopy_iterator
{
	struct all_copy;
	struct all_reference
	{
		Vortex2D& vtx;
		std::pair<int, int>& aflPnl;
		int& num;

		all_reference() = delete;
		all_reference(all_reference const&) = default;
		all_reference(Vortex2D& _vtx, std::pair<int, int>& _aflPnl, int& _num) : vtx(_vtx), aflPnl(_aflPnl), num(_num) {};
		
		~all_reference() = default;

		all_reference(all_reference&& other) = default;

		all_reference& operator=(all_reference&& other)
		{
			vtx = std::move(other.vtx);
			aflPnl = std::move(other.aflPnl);
			num = std::move(other.num);

			return *this;
		}

		friend void swap(all_reference p0, all_reference p1)
		{
			std::swap(p0.vtx, p1.vtx);
			std::swap(p0.aflPnl, p1.aflPnl);
			std::swap(p0.num, p1.num);
		}

		all_reference& operator= (all_copy const& p);

		all_reference& operator= (all_copy&& p);
	};

	struct all_copy
	{
		Vortex2D vtx;
		std::pair<int, int> aflPnl;
		int num;
		
		all_copy(all_reference const& ptr) : vtx(ptr.vtx), aflPnl(ptr.aflPnl), num(ptr.num) {};
		all_copy(all_reference&& ptr) : vtx(std::move(ptr.vtx)), aflPnl(std::move(ptr.aflPnl)), num(std::move(ptr.num)) {};
		
	
	};

	struct all_iterator : public std::iterator< std::random_access_iterator_tag, all_copy >
	{
	private:
		using ItVtx = std::vector<Vortex2D>::iterator;
		using ItAflPnl = std::vector<std::pair<int, int>>::iterator;
		using ItNum = std::vector<int>::iterator;
		ItVtx iVtx;
		ItAflPnl iAflPnl;
		ItNum iNum;

	public:
		all_iterator(ItVtx _iVtx, ItAflPnl _iAflPnl, ItNum _iNum) : iVtx(_iVtx), iAflPnl(_iAflPnl), iNum(_iNum) {};
		
		all_iterator(all_iterator const&) = default;            
		all_iterator& operator= (all_iterator const&) = default;
		~all_iterator() = default;

		using reference = all_reference;
		using pointer = all_reference;

		void swap(all_iterator& other)                          
		{
			std::swap(iVtx, other.iVtx);
			std::swap(iAflPnl, other.iAflPnl);
			std::swap(iNum, other.iNum);
		}

		all_reference operator* ()
		{
			return { *iVtx, *iAflPnl, *iNum };
		}

		all_reference const operator* () const 
		{
			return { *iVtx, *iAflPnl, *iNum };
		}

		all_reference operator-> ()
		{
			return { *iVtx, *iAflPnl, *iNum };
		}

		all_reference const operator-> () const
		{
			return { *iVtx, *iAflPnl, *iNum };
		}

		bool operator==(all_iterator const& other) const        
		{
			return iVtx  == other.iVtx;  
		}

		bool operator!=(all_iterator const& other) const        
		{
			return iVtx != other.iVtx; 
		}

		all_iterator& operator++()
		{
			++iVtx;
			++iAflPnl;
			++iNum;
			return *this;
		}

		all_iterator operator++(int)                            // *++r
		{
			all_iterator temp(*this);
			++(*this);
			return temp;
		}

		all_iterator() = default;                            

		all_iterator& operator--()                              // --r
		{
			--iVtx;
			--iAflPnl;
			--iNum;
			return *this;
		}
		all_iterator operator--(int)                            // r--
		{
			all_iterator temp(*this);
			--(*this);
			return temp;
		}
		
		all_iterator& operator+=(difference_type p)             // r += n
		{
			iVtx += p;
			iAflPnl += p;
			iNum += p;
			return *this;
		}

		all_iterator operator+(difference_type p) const         // a + n
		{
			all_iterator temp(*this);
			temp += p;
			return temp;
		}
		
		friend all_iterator operator+(difference_type p, all_iterator temp)         // n + a
		{
			temp += p;
			return temp;
		}

		all_iterator& operator-=(difference_type p)             // r -= n
		{
			iVtx -= p;
			iAflPnl -= p;
			iNum -= p;
			return *this;
		}
		all_iterator operator-(difference_type p) const         // a - n
		{
			all_iterator temp(*this);
			temp -= p;
			return temp;
		}

		difference_type operator-(all_iterator const& p) const       // b - a
		{
			return iVtx - p.iVtx;  
		}

		all_reference operator[](difference_type p)             // a[n]
		{
			return *(*this + p);
		}
		all_reference const operator[](difference_type p) const // a[n]
		{
			return *(*this + p);
		}

		bool operator<(all_iterator const& p) const             // a < b
		{
			return iVtx < p.iVtx;   
		}
		bool operator>(all_iterator const& p) const             // a > b
		{
			return iVtx > p.iVtx;   
		}
		bool operator>=(all_iterator const& p) const            // a >= b
		{
			return iVtx >= p.iVtx; 
		}
		bool operator<=(all_iterator const& p) const            // a >= b
		{
			return iVtx <= p.iVtx;  
		}

		const Vortex2D& getVtx() const
		{
			return *iVtx;
		}

		int getNum() const
		{
			return *iNum;
		}

		const std::pair<int, int>& getAflPnl() const
		{
			return *iAflPnl;
		}

	};

	bool LessXrr(all_reference const& lhs, all_reference const& rhs);
	bool LessYrr(all_reference const& lhs, all_reference const& rhs);

	bool LessXcr(all_copy const& lhs, all_reference const& rhs);
	bool LessYcr(all_copy const& lhs, all_reference const& rhs);

	bool LessXrc(all_reference const& lhs, all_copy const& rhs);
	bool LessYrc(all_reference const& lhs, all_copy const& rhs);

	bool LessXcc(all_copy const& lhs, all_copy const& rhs);
	bool LessYcc(all_copy const& lhs, all_copy const& rhs);
}

using VM2D::PointsCopy;
#endif