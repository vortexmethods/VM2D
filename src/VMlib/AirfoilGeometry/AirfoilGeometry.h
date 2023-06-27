#ifndef AIRFOILGEOMETRY_H_
#define AIRFOILGEOMETRY_H_

#include "numvector.h"
#include "Point2D.h"

namespace VMlib
{
	struct GeomPoint
	{
	//private:
		Point2D r;
		std::string type;
	//public:

		GeomPoint() {};

		GeomPoint(const Point2D& _r, const std::string& _type)
			: r(_r), type(_type) { }

		~GeomPoint() {};
	};
}

using VMlib::GeomPoint;

#endif