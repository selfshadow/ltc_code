struct mat33
{
    operator glm::mat3() const
    {
        return glm::mat3(m[0], m[1], m[2], m[3], m[4], m[5], m[6], m[7], m[8]);
    }

    double m[9];
};


#include "ltc.inc"

#if 1
mat3 M_GGX(const float cosTheta, const float alpha)
{
	int t = std::max<int>(0, std::min<int>(size-1, (int)floorf(cosTheta * size)));
	int a = std::max<int>(0, std::min<int>(size-1, (int)floorf(sqrtf(alpha) * size)));

	return tabM[a + t*size];
}

float magnitude_GGX(const float cosTheta, const float alpha)
{
	int t = std::max<int>(0, std::min<int>(size-1, (int)floorf(cosTheta * size)));
	int a = std::max<int>(0, std::min<int>(size-1, (int)floorf(sqrtf(alpha) * size)));

	return tabMagnitude[a + t*size];
}
#endif
