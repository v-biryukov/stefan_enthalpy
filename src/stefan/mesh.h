#pragma once

#include <math.h>
#include <functional>
#include "../math/Spaces.h"


template <typename Space>
class Mesh
{
	SPACE_TYPEDEFS;

	IndexVector stride;
	AABB aabb;

	Vector h;

	using MaterialParamsGetter = std::function<material_info<Scalar>(const Vector&)>;
	MaterialParamsGetter materialParamsGetter;

public:

	Mesh() = default;

	Mesh(const IndexVector& stride, const AABB& aabb,
		MaterialParamsGetter materialParamsGetter)
		: stride(stride), aabb(aabb), materialParamsGetter(materialParamsGetter)
	{
		auto size = aabb.Size();
		h = size.ComponentDivide(stride - Vector::one());
	}

	material_info<Scalar> get_material(const IndexVector& cellIndex) const
	{
		return materialParamsGetter(aabb.boxPoint1 + h.ComponentMultiply(cellIndex));
	}

	IndexVector getStride() const { return stride; }
	inline Vector get_h() const { return h; }
	inline IndexType get_number_of_nodes() const { return stride.GetVolume(); }
};