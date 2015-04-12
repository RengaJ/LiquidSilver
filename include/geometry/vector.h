#ifndef GEOMETRY_H
#define GEOMETRY_H

#include <math.h>
#include <assert.h>

// The LS_Vector class represents a direction in the 3D-coordinate space. This
// class should not be used to represent a point (see LS_Point).
class LS_Vector
{
	public:
		// ==========================================================
		// Public Data
		// ==========================================================

		// X-Coordinate
		float x;

		// Y-Coordinate
		float y;

		// Z-Coordinate
		float z;

		// ==========================================================
		// Public Methods
		// ==========================================================

		// Default Constructor
		// ---------------------------------------------
		// Creates and initializes an instance of a LS_Vector
		// to have no direction.
		// ---------------------------------------------
		LS_Vector()
		{
			x = 0.0f;
			y = 0.0f;
			z = 0.0f;
		}

		// Full Constructor
		// ---------------------------------------------
		// Creates and initializes an instance of a
		// LS_Vector pointing in the given (x,y,z) direction
		// ---------------------------------------------
		LS_Vector(const float _x, const float _y, const float _z)
		{
			x = _x;
			y = _y;
			z = _z;

			// Ensure that all of the given values are valid
			// (they are all numbers)
			assert(!hasNaNs());
		}

		// Copy Constructor
		// ---------------------------------------------
		// Creates and initializes an instance of a
		// LS_Vector that points in the same direction as
		// another LS_Vector
		// ---------------------------------------------
		LS_Vector(LS_Vector& direction)
		{
			x = direction.x;
			y = direction.y;
			z = direction.z;

			// Ensure that all of the given values are valid
			// (they are all numbers)
			assert(!hasNaNs);
		}

		// ---------------------------------------------
		// Determines if the current vector has any value
		// that is considered to be "not a number (nan)".
		// ---------------------------------------------
		bool hasNaNs() const
		{
			return (isnan(x) || isnan(y) || isnan(z));
		}

		// Addition Operator
		// ---------------------------------------------
		// Add the current LS_Vector and another LS_Vector
		// together to create a third LS_Vector that is the
		// sum of both vectors.
		// ---------------------------------------------
		LS_Vector operator+(const LS_Vector &v) const
		{
			return LS_Vector(x + v.x, y + v.y, z + v.z);
		}

		// Compound Addition Operator
		// ---------------------------------------------
		// Add the given LS_Vector to the current LS_Vector
		// and return the current LS_Vector
		// ---------------------------------------------
		LS_Vector& operator+=(const LS_Vector &v)
		{
			x += v.x;
			y += v.y;
			z += v.z;

			return *this;
		}

		// Subtraction Operator
		// ---------------------------------------------
		// Subtract the current LS_Vector and the given LS_Vector
		// from each other to create a third LS_Vector that is the
		// difference between the two vectors
		// ---------------------------------------------
		LS_Vector operator-(const LS_Vector &v) const
		{
			return LS_Vector(x - v.x, y - v.y, z - v.z);
		}

		// Compound Subtraction Operator
		// ---------------------------------------------
		// Subtract the given LS_Vector from the current
		// LS_Vector and return the current LS_Vector
		// ---------------------------------------------
		LS_Vector& operator-=(const LS_Vector &v)
		{
			x -= v.x;
			y -= v.y;
			z -= v.z;

			return *this;
		}

		// Scalar Multiplication
		// ---------------------------------------------
		// Multiply the current LS_Vector by a given scalar
		// and create a new LS_Vector that is the scaled result
		// ---------------------------------------------
		LS_Vector operator*(float scalar)
		{
			return LS_Vector(scalar * x, scalar * y, scalar * z);
		}

		// Scalar Division
		// ---------------------------------------------
		// Divide the current LS_Vector by a given scalar
		// and create a new LS_Vector that is the scaled result
		// ---------------------------------------------
		LS_Vector operator/(float scalar) const
		{
			assert(scalar != 0);

			float reciprocal = 1.f / scalar;
			// Note how this operator uses multiplication and only one division...
			// This is due to the fact that mutliplication is faster than division,
			// and optimizes the division operation
			return LS_Vector(x * reciprocal, y * reciprocal, z * reciprocal);
		}
		
		// Compound Scalar Division
		// ---------------------------------------------
		// Divide the current LS_Vector by the given scalar
		// and return the current LS_Vector
		// ---------------------------------------------
		LS_Vector& operator/=(float scalar)
		{
			assert(scalar != 0);

			float reciprocal = 1.f / scalar;
			// Note how this operator uses multiplication and only one division...
			// This is due to the fact that multiplication is faster than division,
			// and optimizes the division operation
			x *= reciprocal;
			y *= reciprocal;
			z *= reciprocal;

			return *this;
		}

		// Unary Negation
		// ---------------------------------------------
		// Inverts the direction of the LS_Vector 180 degrees (PI Radians)
		// and returns a LS_Vector pointing in this new direction
		// ---------------------------------------------
		LS_Vector operator-() const
		{
			return LS_Vector(-x, -y, -z);
		}

		// Index Notation (I)
		// ---------------------------------------------
		// Gain introspection into the LS_Vector component by
		// using standard array-index notation.
		// ---------------------------------------------
		float operator[](int i) const
		{
			assert((i >= 0) && (i <= 2));

			return (&x)[i];
		}

		// Index Notation (II)
		// ---------------------------------------------
		// Gain introspection into the LS_Vector component by
		// using standard array-index notation.
		// ---------------------------------------------
		float &operator[](int i)
		{
			assert((i >= 0) && (i <= 2));

			return (&x)[i];
		}

		// Length (Squared)
		// ---------------------------------------------
		// Calculate and return the squared length of the
		// LS_Vector. This is faster than calculating the
		// actual length of the LS_Vector
		// ---------------------------------------------
		float lengthSquared() const
		{
			return (x * x) + (y * y) + (z * z);
		}

		// Length
		// ---------------------------------------------
		// Calculate and return the length of the LS_Vector
		// ---------------------------------------------
		float length() const
		{
			return sqrtf(lengthSquared());
		}
};

// ---------------------------------------------
// |        LS_Vector Inline Functions         |
// ---------------------------------------------

// LS_Vector: Scalar Multiplication
// ---------------------------------------------
// Multiply the given LS_Vector by the given scalar
// and create a new LS_Vector that is the scaled result
// ---------------------------------------------
inline LS_Vector operator*(float scalar, const LS_Vector &v)
{
	return LS_Vector(scalar * v.x, scalar * v.y, scalar * v.z);
}

// LS_Vector: Dot Product
// ---------------------------------------------
// Perform a Dot Operation on two LS_Vectors v and w
// such that:
//    S = (v_x * w_x) + (v_y * w_y) + (v_z * w_z)
// ---------------------------------------------
inline float Dot(const LS_Vector &v1, const LS_Vector &v2)
{
	return (v1.x * v2.x) + (v1.y * v2.y) + (v1.z * v2.z);
}

// LS_Vector: Absolute-Dot Product
// ---------------------------------------------
// Perform a Dot Operation on two LS_Vectors v and w
// and obtain the absolute value of the result:
//    S = |(v_x * w_x) + (v_y * w_y) + (v_z * w_z)|
// ---------------------------------------------
inline float AbsDot(const LS_Vector &v1, const LS_Vector &v2)
{
	return fabsf(Dot(v1, v2));
}

// LS_Vector: Cross Product
// ---------------------------------------------
// Perform a Cross Production operation on two
// LS_Vectors v and w such that:
//  V = { (v_y * w_z) - (v_z * w_y),
//        (v_z * w_x) - (v_x * w_z),
//        (v_x * w_y) - (v_y * w_x) }
// ---------------------------------------------
inline LS_Vector Cross(const LS_Vector &v1, const LS_Vector &v2)
{
	return LS_Vector(
		(v1.y * v2.z) - (v1.z * v2.y),
		(v1.z * v2.x) - (v1.x * v2.z),
		(v1.x * v2.y) - (v1.y * v2.x));
}

// LS_Vector: Normalize
// ---------------------------------------------
// Normalize a LS_Vector such that its length is
// equal to 1
// ---------------------------------------------
inline LS_Vector Normalize(const LS_Vector &v)
{
	return (v / v.length());
}

// LS_Vector: Coordinate System
// ---------------------------------------------
// Calculate and set the coordinate system for three LS_Vectors,
// given a single LS_Vector
// ---------------------------------------------
inline void CoordinateSystem(const LS_Vector &v1, LS_Vector* v2, LS_Vector *v3)
{
	if (fabsf(v1.x) > fabs(v1.y))
	{
		float reciprocal = 1.f / sqrt((v1.x * v1.x) + (v1.z * v1.z));
		*v2 = LS_Vector(-v1.z * reciprocal, 0.f, v1.x * reciprocal);
	}
	else
	{
		float reciprocal = 1.f / sqrt((v1.y * v1.y) + (v1.z * v1.z));
		*v2 = LS_Vector(0.f, v1.z * reciprocal, -v1.y * reciprocal);
	}

	*v3 = Cross(v1, *v2);
}
#endif