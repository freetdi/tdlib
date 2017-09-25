
#ifndef HAVE_BOOL_H
#define HAVE_BOOL_H
class BOOL{ //
public:
	BOOL() : value_(bool())
	{
	}
	/* explicit */ BOOL(bool const& t): value_(t) {}
	// /* explicit */ operator bool&() { return value_; }
	/* explicit */ operator bool() const { return value_; }
private:
	char value_;
};
#endif
