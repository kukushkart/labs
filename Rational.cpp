#include "Rational.h"

int gcd(int a, int b) {
	while (b != 0) {
		int temp = b;
		b = a % b;
		a = temp;
	}
	return a;
}
void Rational::redaction()
{
	int gcd_ = gcd(numerator, denominator);
	numerator /= gcd_;
	denominator /= gcd_;
	if (denominator < 0) {
		numerator *= -1;
		denominator *= -1;
	}
}



Rational Rational::operator+(const Rational& rhs)
{
	return Rational(numerator * rhs.denominator + denominator * rhs.numerator, denominator * rhs.denominator);
}



Rational Rational::operator-(const Rational& rhs)
{
	return Rational(this->numerator * rhs.denominator - this->denominator * rhs.numerator, this->denominator * rhs.denominator);
}


Rational Rational::operator*(const Rational& rhs) const
{
	return Rational(numerator * rhs.numerator, denominator * rhs.denominator);
}




Rational Rational::operator/(const Rational& rhs)
{
	return Rational(numerator * rhs.denominator, denominator * rhs.numerator);
}

Rational Rational::operator/(int n)
{
	return Rational(numerator, denominator * n);
}



std::istream& operator>>(std::istream& in, Rational& a)
{
	in >> a.numerator;
	in.ignore('/');
	in >> a.denominator;
	return in;
}

std::ostream& operator<<(std::ostream& os, const Rational& a)
{
	os << a.numerator << " / " << a.denominator << std::endl;
	return os;
}

Rational abs(const Rational& r) {
	return Rational(std::abs(r.numerator), std::abs(r.denominator));
}


Rational::Rational(int new_numerator, int new_denominator) :numerator(new_numerator), denominator(new_denominator) {
	if (denominator == 0) {
		throw "Dominator cannot be zero.";
	}
	(*this).redaction();
}



Rational::Rational(int num) :numerator(num), denominator(1) {}

Rational& Rational::operator*=(const Rational& rhs)
{
	(*this).numerator = (*this).numerator * rhs.numerator;
	(*this).denominator = (*this).denominator * rhs.denominator;
	return *this;
}

Rational& Rational::operator*=(int n)
{
	(*this).numerator = (*this).numerator * n;
	return *this;
}

Rational& Rational::operator=(const Rational& rhs)
{
	if (this == &rhs) return *this;
	numerator = rhs.numerator;
	denominator = rhs.denominator;
	return *this;

}

Rational& Rational::operator=(const int n)
{
	numerator = n;
	denominator = 1;
	return *this;
}

Rational& Rational::operator+=(const Rational& rhs)
{
	(*this).denominator = (*this).denominator * rhs.denominator;
	(*this).numerator = (*this).numerator * rhs.denominator + (*this).denominator * rhs.numerator;
	(*this).redaction();
	return *this;
}
Rational& Rational::operator-=(const Rational& rhs)
{
	(*this).denominator = (*this).denominator * rhs.denominator;
	(*this).numerator = (*this).numerator * rhs.denominator - (*this).denominator * rhs.numerator;
	(*this).redaction();
	return *this;
}

Rational& Rational::operator/=(const Rational& rhs)
{
	(*this).numerator *= rhs.denominator;
	(*this).denominator *= rhs.numerator;
	return *this;
}

bool Rational::operator==(const Rational& rhs) const
{
	if (rhs.numerator == (*this).numerator && rhs.denominator == (*this).denominator) {
		return true;
	}
	else {
		return false;
	}
}

bool Rational::operator>(const Rational& rhs) const
{
	if (rhs.denominator != (*this).denominator) {
		if ((*this).numerator * rhs.denominator > rhs.numerator * (*this).denominator) {
			return true;
		}
		else {
			return false;
		}
	}
	else {
		if ((*this).numerator > rhs.denominator) {
			return true;
		}
		else { return false; }
	}
}

bool Rational::operator!=(const Rational& rhs) const
{
	if (!(rhs == (*this))) {
		return true;
	}
	return false;
}

