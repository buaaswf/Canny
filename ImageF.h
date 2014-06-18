#pragma once
class ImageF
{
public:
	ImageF(void);
	~ImageF(void);
	friend ImageF operator/(ImageF x,float s);
};

