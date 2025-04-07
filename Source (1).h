#pragma once
#ifndef RATIONAL_H
#define RATIONAL_H
#include <fstream>
class Rational
{
private:
	int numerator;
	int denominator;
public:
	void redaction();
	Rational() :numerator(0), denominator(1) {}
	Rational(int new_numerator, int new_denominator);
	Rational operator+(const Rational& rhs);
	Rational operator-(const Rational& rhs);
	Rational operator*(const Rational& rhs) const;
	Rational operator/(const Rational& rhs);
	Rational operator/(int n);
	Rational(int num);
	Rational& operator*=(const Rational& rhs);
	Rational& operator*=(int n);
	Rational& operator=(const Rational& rhs);
	Rational& operator=(const int n);
	Rational& operator+=(const Rational& rhs);
	Rational& operator-=(const Rational& rhs);
	Rational& operator/=(const Rational& rhs);
	bool operator==(const Rational& rhs) const;
	bool operator>(const Rational& rhs) const;
	bool operator!=(const Rational& rhs) const;
	friend std::istream& operator>>(std::istream& in, Rational& a);
	friend std::ostream& operator<<(std::ostream& os, const Rational& a);
	friend Rational abs(const Rational& r);
};

#endif


