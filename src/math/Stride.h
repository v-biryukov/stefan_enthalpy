#pragma once

#include "Spaces.h"

template <typename Space>
class AllNodes;

template <typename Space>
class StrideIterator 
{
public:
	SPACE_TYPEDEFS;
	using Stride = IndexVector;

	StrideIterator(const Stride& value, 
		const Stride& maxValue):
		value(value),
		maxValue(maxValue)
	{
	}

	StrideIterator& operator++()
	{
		++value[0];
		for (int i = 0; i < Space::Dimension; ++i)
		{
			if (value[i] == maxValue[i])
			{
				if (i + 1 < Space::Dimension)
				{
					value[i] = 0;
					value[i + 1]++;
				}
			}
		}
		return *this;
	}

	bool operator!=(const StrideIterator& other) const
	{
		for (int i = 0; i < Space::Dimension; ++i)
			if (value[i] != other.value[i])
				return true;
		return false;
	}

	Stride operator*() const
	{
		return value;
	}

	friend class AllNodes<Space>;

private:
	Stride value;
	Stride maxValue;
};

template <typename Space>
class AllNodes
{
public:
	SPACE_TYPEDEFS;
	using Stride = IndexVector;

	AllNodes(const Stride& boundary) : boundary(boundary)
	{
	}

	StrideIterator<Space> begin() const
	{
		return { Stride::zero(), boundary };
	}

	StrideIterator<Space> end() const
	{
		Stride e = Stride::zero();
		e[Space::Dimension - 1] = boundary[Space::Dimension - 1];
		return { e, boundary };
	}

private:
	Stride boundary;
};
