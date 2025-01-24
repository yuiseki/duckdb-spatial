#pragma once

#include <cstdint>
#include <cstring>
#include <string>
#include <cmath>

// Assert macro
#ifndef SGL_ASSERT
#ifdef NDEBUG
#define SGL_ASSERT(x) ((void)0)
#else
#include <cassert>
#define SGL_ASSERT(x) assert(x)
#endif
#endif

namespace sgl {

struct allocator {
	virtual void *alloc(size_t size) = 0;
	virtual void dealloc(void *ptr, size_t size) = 0;
	virtual void *realloc(void *ptr, size_t old_size, size_t new_size) = 0;
	virtual ~allocator() = default;
};

struct vertex_xy {
	double x;
	double y;

	bool operator==(const vertex_xy &other) const {
		return x == other.x && y == other.y;
	}
};

struct vertex_xyzm {
	double x;
	double y;
	double zm;
	double m;

	bool operator==(const vertex_xyzm &other) const {
		return x == other.x && y == other.y && zm == other.zm && m == other.m;
	}
};

struct box_xy {
	vertex_xy min;
	vertex_xy max;

	static box_xy smallest() {
		constexpr auto dmax = std::numeric_limits<double>::max();
		constexpr auto dmin = std::numeric_limits<double>::lowest();
		return {
		    {dmax, dmax},
		    {dmin, dmin},
		};
	}
};

struct box_xyzm {
	vertex_xyzm min;
	vertex_xyzm max;

	static box_xyzm smallest() {
		constexpr auto dmax = std::numeric_limits<double>::max();
		constexpr auto dmin = std::numeric_limits<double>::lowest();
		return {
		    {dmax, dmax, dmax, dmax},
		    {dmin, dmin, dmin, dmin},
		};
	}
};

enum class geometry_type : uint8_t {
	INVALID = 0,
	POINT,
	LINESTRING,
	POLYGON,
	MULTI_POINT,
	MULTI_LINESTRING,
	MULTI_POLYGON,
	MULTI_GEOMETRY,
};

enum class vertex_type : uint8_t {
	XY = 0,
	XYZ = 1,
	XYM = 2,
	XYZM = 3,
};

class geometry {
private:
	// clang-format off
	geometry_type	type = geometry_type::INVALID;
	uint8_t			flag = 0;
	uint16_t		padd = 0;
	uint32_t		size = 0;
	void*			data = nullptr;
	geometry*		next = nullptr;
	geometry*		prnt = nullptr;
	// clang-format on
public:
	geometry() = default;
	explicit geometry(const geometry_type type, const bool has_z = false, const bool has_m = false) : type(type) {
		set_z(has_z);
		set_m(has_m);
	}

	geometry_type get_type() const;
	void set_type(geometry_type type);

	bool has_z() const;
	bool has_m() const;
	bool set_z(bool value);
	bool set_m(bool value);

	bool is_single_part() const;
	bool is_multi_part() const;
	bool is_collection() const;

	uint32_t get_count() const;
	void set_count(uint32_t count);
	bool is_empty() const;

	const geometry *get_last_part() const;
	const geometry *get_first_part() const;
	const geometry *get_nth_part(uint32_t n) const;
	const geometry *get_next() const;
	const geometry *get_parent() const;

	geometry *get_last_part();
	geometry *get_first_part();
	geometry *get_nth_part(uint32_t n);
	geometry *get_next();
	geometry *get_parent();

	void append_part(geometry *part);

	const uint8_t *get_vertex_data() const;
	uint8_t *get_vertex_data();
	void set_vertex_data(const uint8_t *data, uint32_t size);
	void set_vertex_data(const char *data, uint32_t size);

	size_t get_vertex_size() const;
	vertex_xy get_vertex_xy(uint32_t n) const;
	vertex_xyzm get_vertex_xyzm(uint32_t n) const;

	static std::string type_to_string(geometry_type type);
};

} // namespace sgl

//--------------------------------------------------------------------------
// Implementation
//--------------------------------------------------------------------------

namespace sgl {

inline geometry_type geometry::get_type() const {
	return type;
}

inline void geometry::set_type(const geometry_type type) {
	this->type = type;
}

inline bool geometry::has_z() const {
	return flag & 0x01;
}

inline bool geometry::has_m() const {
	return flag & 0x02;
}

inline bool geometry::set_z(const bool value) {
	if (value) {
		flag |= 0x01;
	} else {
		flag &= ~0x01;
	}
	return value;
}

inline bool geometry::set_m(const bool value) {
	if (value) {
		flag |= 0x02;
	} else {
		flag &= ~0x02;
	}
	return value;
}

inline bool geometry::is_single_part() const {
	return type == geometry_type::POINT || type == geometry_type::LINESTRING;
}

inline bool geometry::is_multi_part() const {
	return type >= geometry_type::POLYGON && type <= geometry_type::MULTI_GEOMETRY;
}

inline bool geometry::is_collection() const {
	return type >= geometry_type::MULTI_POINT && type <= geometry_type::MULTI_GEOMETRY;
}

inline uint32_t geometry::get_count() const {
	return size;
}

inline void geometry::set_count(const uint32_t count) {
	size = count;
}

inline bool geometry::is_empty() const {
	return size == 0;
}

inline const geometry *geometry::get_last_part() const {
	SGL_ASSERT(is_multi_part() || type == geometry_type::INVALID);
	const auto tail = static_cast<geometry *>(data);
	return tail;
}

inline const geometry *geometry::get_first_part() const {
	const auto tail = get_last_part();
	return tail ? tail->next : nullptr;
}

inline const geometry *geometry::get_nth_part(uint32_t n) const {
	if (size == 0) {
		SGL_ASSERT(data == nullptr);
		return nullptr;
	}

	auto part = get_first_part();
	SGL_ASSERT(part != nullptr);

	for (uint32_t i = 0; i < n; i++) {
		part = part->next;
		SGL_ASSERT(part != nullptr);
	}

	return part;
}

inline const geometry *geometry::get_next() const {
	return next;
}

inline const geometry *geometry::get_parent() const {
	return prnt;
}

inline geometry *geometry::get_last_part() {
	return const_cast<geometry *>(static_cast<const geometry *>(this)->get_last_part());
}
inline geometry *geometry::get_first_part() {
	return const_cast<geometry *>(static_cast<const geometry *>(this)->get_first_part());
}
inline geometry *geometry::get_nth_part(uint32_t n) {
	return const_cast<geometry *>(static_cast<const geometry *>(this)->get_nth_part(n));
}
inline geometry *geometry::get_next() {
	return const_cast<geometry *>(static_cast<const geometry *>(this)->get_next());
}
inline geometry *geometry::get_parent() {
	return const_cast<geometry *>(static_cast<const geometry *>(this)->get_parent());
}

inline void geometry::append_part(geometry *part) {
	SGL_ASSERT(is_multi_part() || type == geometry_type::INVALID);
	SGL_ASSERT(part != nullptr);

	const auto tail = static_cast<geometry *>(data);

	if (tail == nullptr) {
		SGL_ASSERT(size == 0);
		part->next = part;
	} else {
		SGL_ASSERT(size != 0);
		const auto head = tail->next;
		tail->next = part;
		part->next = head;
	}

	part->prnt = this;
	data = part;
	size++;
}

inline const uint8_t *geometry::get_vertex_data() const {
	SGL_ASSERT(is_single_part() || type == geometry_type::INVALID);
	return static_cast<const uint8_t *>(data);
}

inline uint8_t *geometry::get_vertex_data() {
	SGL_ASSERT(is_single_part() || type == geometry_type::INVALID);
	return static_cast<uint8_t *>(data);
}

inline void geometry::set_vertex_data(const uint8_t *data, uint32_t size) {
	SGL_ASSERT(is_single_part() || type == geometry_type::INVALID);
	// Points can have at most one vertex
	SGL_ASSERT(type != geometry_type::POINT || size < 2);
	this->data = const_cast<uint8_t *>(data);
	this->size = size;
}

inline void geometry::set_vertex_data(const char *data, uint32_t size) {
	set_vertex_data(reinterpret_cast<const uint8_t *>(data), size);
}

inline size_t geometry::get_vertex_size() const {
	return sizeof(double) * (2 + has_z() + has_m());
}

inline vertex_xy geometry::get_vertex_xy(const uint32_t n) const {
	SGL_ASSERT(is_single_part() || type == geometry_type::INVALID);
	SGL_ASSERT(n < size);

	const auto vertex_stride = get_vertex_size();
	const auto vertex_buffer = get_vertex_data();
	const auto vertex_offset = vertex_buffer + vertex_stride * n;

	vertex_xy vertex = {0};
	memcpy(&vertex, vertex_offset, sizeof(vertex_xy));
	return vertex;
}

inline vertex_xyzm geometry::get_vertex_xyzm(const uint32_t n) const {
	SGL_ASSERT(is_single_part() || type == geometry_type::INVALID);
	SGL_ASSERT(n < size);

	const auto vertex_stride = get_vertex_size();
	const auto vertex_buffer = get_vertex_data();
	const auto vertex_offset = vertex_buffer + vertex_stride * n;

	vertex_xyzm vertex = {0};
	memcpy(&vertex, vertex_offset, vertex_stride);
	return vertex;
}

inline std::string geometry::type_to_string(const geometry_type type) {
	switch (type) {
	case geometry_type::POINT:
		return "POINT";
	case geometry_type::LINESTRING:
		return "LINESTRING";
	case geometry_type::POLYGON:
		return "POLYGON";
	case geometry_type::MULTI_POINT:
		return "MULTIPOINT";
	case geometry_type::MULTI_LINESTRING:
		return "MULTILINESTRING";
	case geometry_type::MULTI_POLYGON:
		return "MULTIPOLYGON";
	case geometry_type::MULTI_GEOMETRY:
		return "GEOMETRYCOLLECTION";
	default:
		return "INVALID";
	}
}

} // namespace sgl

//--------------------------------------------------------------------------
// Operations
//--------------------------------------------------------------------------

namespace sgl {

namespace point {
inline geometry make_empty(bool has_z = false, bool has_m = false) {
	return geometry(geometry_type::POINT, has_z, has_m);
}
} // namespace point

namespace linestring {
inline geometry make_empty(bool has_z = false, bool has_m = false) {
	return geometry(geometry_type::LINESTRING, has_z, has_m);
}

inline bool is_closed(const geometry *geom) {
	SGL_ASSERT(geom->get_type() == geometry_type::LINESTRING);

	if (geom->get_count() < 2) {
		return false;
	}

	const auto first = geom->get_vertex_xyzm(0);
	const auto last = geom->get_vertex_xyzm(geom->get_count() - 1);
	return first == last;
}

inline double signed_area(const geometry *geom) {
	SGL_ASSERT(geom->get_type() == geometry_type::LINESTRING);
	SGL_ASSERT(is_closed(geom));

	const auto count = geom->get_count();

	if (count < 3) {
		return 0.0;
	}

	const auto vertex_data = geom->get_vertex_data();
	const auto vertex_size = geom->get_vertex_size();

	auto area = 0.0;

	double x0 = 0.0;
	double x1 = 0.0;
	double y1 = 0.0;
	double y2 = 0.0;

	const auto x_data = vertex_data;
	const auto y_data = vertex_data + sizeof(double);

	memcpy(&x0, x_data, sizeof(double));

	for (uint32_t i = 1; i < count - 1; i++) {
		memcpy(&x1, x_data + (i + 0) * vertex_size, sizeof(double));
		memcpy(&y1, y_data + (i + 1) * vertex_size, sizeof(double));
		memcpy(&y2, y_data + (i - 1) * vertex_size, sizeof(double));

		area += (x1 - x0) * (y2 - y1);
	}

	return area * 0.5;
}

inline double length(const geometry *geom) {
	SGL_ASSERT(geom->get_type() == geometry_type::LINESTRING);

	const auto count = geom->get_count();
	if (count < 2) {
		return 0.0;
	}

	const auto vertex_data = geom->get_vertex_data();
	const auto vertex_size = geom->get_vertex_size();

	auto length = 0.0;

	vertex_xy prev = {0};
	vertex_xy next = {0};

	memcpy(&prev, vertex_data, sizeof(vertex_xy));

	for (uint32_t i = 1; i < count; i++) {
		memcpy(&next, vertex_data + i * vertex_size, sizeof(vertex_xy));
		const auto dx = next.x - prev.x;
		const auto dy = next.y - prev.y;
		length += std::hypot(dx, dy);
		prev = next;
	}

	return length;
}

} // namespace linestring

namespace polygon {

inline geometry make_empty(bool has_z = false, bool has_m = false) {
	return geometry(geometry_type::POLYGON, has_z, has_m);
}

inline double area(const geometry *geom) {
	SGL_ASSERT(geom->get_type() == geometry_type::POLYGON);
	double area = 0.0;

	auto part = geom->get_first_part();
	if (!part) {
		return area;
	}

	area += std::abs(linestring::signed_area(part));

	while (part != geom->get_last_part()) {
		part = part->get_next();
		area -= std::abs(linestring::signed_area(part));
	}

	return area;
}

inline double perimeter(const geometry *geom) {
	SGL_ASSERT(geom->get_type() == geometry_type::POLYGON);

	const auto tail = geom->get_last_part();
	if (!tail) {
		return 0.0;
	}

	double perimeter = 0.0;
	auto part = tail;
	do {
		part = part->get_next();
		perimeter += linestring::length(part);
	} while (part != tail);

	return perimeter;
}

} // namespace polygon

namespace multi_point {

inline geometry make_empty(bool has_z = false, bool has_m = false) {
	return geometry(geometry_type::MULTI_POINT, has_z, has_m);
}

}; // namespace multi_point

namespace multi_linestring {
inline geometry make_empty(bool has_z = false, bool has_m = false) {
	return geometry(geometry_type::MULTI_LINESTRING, has_z, has_m);
}

inline bool is_closed(const geometry *geom) {
	SGL_ASSERT(geom->get_type() == geometry_type::MULTI_LINESTRING);

	const auto tail = geom->get_last_part();
	if (!tail) {
		return false;
	}

	auto part = tail;
	do {
		part = part->get_next();
		if (!linestring::is_closed(part)) {
			return false;
		}
	} while (part != tail);

	return true;
}

inline double length(const geometry *geom) {
	SGL_ASSERT(geom->get_type() == geometry_type::MULTI_LINESTRING);

	const auto tail = geom->get_last_part();
	if (!tail) {
		return 0.0;
	}

	double length = 0.0;
	auto part = tail;
	do {
		part = part->get_next();
		length += linestring::length(part);
	} while (part != tail);

	return length;
}

} // namespace multi_linestring

namespace multi_polygon {
inline geometry make_empty(bool has_z = false, bool has_m = false) {
	return geometry(geometry_type::MULTI_POLYGON, has_z, has_m);
}

inline double area(const geometry *geom) {
	SGL_ASSERT(geom->get_type() == geometry_type::MULTI_POLYGON);

	const auto tail = geom->get_last_part();
	if (!tail) {
		return 0.0;
	}

	double area = 0.0;
	auto part = tail;
	do {
		part = part->get_next();
		area += polygon::area(part);
	} while (part != tail);

	return area;
}

inline double perimeter(const geometry *geom) {
	SGL_ASSERT(geom->get_type() == geometry_type::MULTI_POLYGON);

	const auto tail = geom->get_last_part();
	if (!tail) {
		return 0.0;
	}

	double perimeter = 0.0;
	auto part = tail;
	do {
		part = part->get_next();
		perimeter += polygon::perimeter(part);
	} while (part != tail);

	return perimeter;
}

} // namespace multi_polygon

namespace multi_geometry {

inline geometry make_empty(bool has_z = false, bool has_m = false) {
	return geometry(geometry_type::MULTI_GEOMETRY, has_z, has_m);
}
inline double area(const geometry *geom);
inline double length(const geometry *geom);
inline double perimeter(const geometry *geom);

}; // namespace multi_geometry

namespace ops {

double area(const geometry *geom);
double perimeter(const geometry *geom);
double length(const geometry *geom);
size_t vertex_count(const geometry *geom);
int32_t max_surface_dimension(const geometry *geom);

box_xy extent_xy(const geometry *geom);
void force_zm(allocator &alloc, geometry *geom, bool has_z, bool has_m, double default_z, double default_m);

} // namespace ops

} // namespace sgl

//--------------------------------------------------------------------------
// Implementation
//--------------------------------------------------------------------------

namespace sgl {
namespace ops {

inline double area(const geometry *geom) {
	switch (geom->get_type()) {
	case geometry_type::POLYGON: {
		return polygon::area(geom);
	}
	case geometry_type::MULTI_POLYGON: {
		return multi_polygon::area(geom);
	}
	case geometry_type::MULTI_GEOMETRY: {
		return multi_geometry::area(geom);
	}
	default:
		return 0.0;
	}
}

inline double perimeter(const geometry *geom) {
	switch (geom->get_type()) {
	case geometry_type::POLYGON: {
		return polygon::perimeter(geom);
	}
	case geometry_type::MULTI_POLYGON: {
		return multi_polygon::perimeter(geom);
	}
	case geometry_type::MULTI_GEOMETRY: {
		return multi_geometry::perimeter(geom);
	}
	default:
		return 0.0;
	}
}

inline double length(const geometry *geom) {
	switch (geom->get_type()) {
	case geometry_type::LINESTRING: {
		return linestring::length(geom);
	}
	case geometry_type::MULTI_LINESTRING: {
		return multi_linestring::length(geom);
	}
	case geometry_type::MULTI_GEOMETRY: {
		return multi_geometry::length(geom);
	}
	default:
		return 0;
	}
}

inline size_t vertex_count(const geometry *geom) {
	if (geom->is_single_part()) {
		return geom->get_count();
	} else {
		const auto tail = geom->get_last_part();
		if (!tail) {
			return 0;
		}
		auto part = tail;
		size_t count = 0;
		do {
			part = part->get_next();
			count += part->get_count();
		} while (part != tail);
		return count;
	}
}

inline int32_t max_surface_dimension(const geometry *geom) {
	switch (geom->get_type()) {
	case geometry_type::POINT:
	case geometry_type::MULTI_POINT:
		return 0;
	case geometry_type::LINESTRING:
	case geometry_type::MULTI_LINESTRING:
		return 1;
	case geometry_type::POLYGON:
	case geometry_type::MULTI_POLYGON:
		return 2;
	case geometry_type::MULTI_GEOMETRY: {
		int32_t max_dim = 0;
		const auto tail = geom->get_last_part();
		if (!tail) {
			return 0;
		}
		auto part = tail;
		do {
			part = part->get_next();
			max_dim = std::max(max_dim, max_surface_dimension(part));
		} while (part != tail);
		return max_dim;
	}
	default:
		return 0;
	}
}

template <class F>
inline void visit_vertices(const geometry *geom, F callback) {
	switch (geom->get_type()) {
	case geometry_type::POINT:
	case geometry_type::LINESTRING: {
		auto vertex_data = geom->get_vertex_data();
		if (vertex_data == nullptr) {
			return;
		}
		auto vertex_size = geom->get_vertex_size();
		for (uint32_t i = 0; i < geom->get_count(); i++) {
			callback(vertex_data + i * vertex_size);
		}
	}
		return;
	case geometry_type::POLYGON:
	case geometry_type::MULTI_POINT:
	case geometry_type::MULTI_LINESTRING:
	case geometry_type::MULTI_POLYGON:
	case geometry_type::MULTI_GEOMETRY: {
		const auto tail = geom->get_last_part();
		if (!tail) {
			return;
		}
		auto part = tail;
		do {
			part = part->get_next();
			visit_vertices<F>(part, callback);
		} while (part != tail);
	}
		return;
	default:
		return;
	}
}

// Non-recursive
inline bool try_get_extent_xy(const geometry *geom, box_xy *out) {

	auto part = geom;
	if (part == nullptr) {
		return false;
	}

	box_xy result = box_xy::smallest();

	const auto root = part->get_parent();
	bool has_any_vertices = false;

	while (true) {
		switch (part->get_type()) {
		case geometry_type::POINT:
		case geometry_type::LINESTRING: {
			const auto vertex_count = part->get_count();

			has_any_vertices |= vertex_count > 0;

			for (uint32_t i = 0; i < vertex_count; i++) {
				const auto vertex = part->get_vertex_xy(i);
				result.min.x = std::min(result.min.x, vertex.x);
				result.min.y = std::min(result.min.y, vertex.y);
				result.max.x = std::max(result.max.x, vertex.x);
				result.max.y = std::max(result.max.y, vertex.y);
			}
		} break;
		case geometry_type::POLYGON: {
			if (!part->is_empty()) {
				const auto shell = part->get_first_part();
				const auto vertex_count = shell->get_count();
				has_any_vertices |= vertex_count > 0;
				for (uint32_t i = 0; i < vertex_count; i++) {
					const auto vertex = shell->get_vertex_xy(i);
					result.min.x = std::min(result.min.x, vertex.x);
					result.min.y = std::min(result.min.y, vertex.y);
					result.max.x = std::max(result.max.x, vertex.x);
					result.max.y = std::max(result.max.y, vertex.y);
				}
			}
		} break;
		case geometry_type::MULTI_POINT:
		case geometry_type::MULTI_LINESTRING:
		case geometry_type::MULTI_POLYGON:
		case geometry_type::MULTI_GEOMETRY: {
			if (!part->is_empty()) {
				part = part->get_first_part();
				// continue the outer loop here!
				continue;
			}
		} break;
		default:
			SGL_ASSERT(false);
			return false;
		}

		// Now go up/sideways
		while (true) {
			const auto parent = part->get_parent();
			if (parent == root) {
				if (has_any_vertices) {
					*out = result;
					return true;
				}
				return false;
			}

			if (part != parent->get_last_part()) {
				// Go sideways
				part = part->get_next();
				break;
			}

			// Go up
			part = parent;
		}
	}
}

} // namespace ops
} // namespace sgl

namespace sgl {
namespace multi_geometry {

inline double area(const geometry *geom) {
	SGL_ASSERT(geom->get_type() == geometry_type::MULTI_GEOMETRY);

	const auto tail = geom->get_last_part();
	if (!tail) {
		return 0.0;
	}

	double area = 0.0;
	auto part = tail;
	do {
		part = part->get_next();
		area += ops::area(part);
	} while (part != tail);

	return area;
}

inline double length(const geometry *geom) {
	SGL_ASSERT(geom->get_type() == geometry_type::MULTI_GEOMETRY);

	const auto tail = geom->get_last_part();
	if (!tail) {
		return 0.0;
	}
	double length = 0.0;
	auto part = tail;
	do {
		part = part->get_next();
		length += ops::length(part);
	} while (part != tail);
	return length;
}

inline double perimeter(const geometry *geom) {
	SGL_ASSERT(geom->get_type() == geometry_type::MULTI_GEOMETRY);

	const auto tail = geom->get_last_part();
	if (!tail) {
		return 0.0;
	}
	double perimeter = 0.0;
	auto part = tail;
	do {
		part = part->get_next();
		perimeter += ops::perimeter(part);
	} while (part != tail);
	return perimeter;
}

} // namespace multi_geometry

namespace util {

inline double haversine_distance(const double lat1_p, const double lon1_p, const double lat2_p, const double lon2_p) {
	// Radius of the earth in km
	constexpr auto R = 6371000.0;
	constexpr auto PI = 3.14159265358979323846;

	// Convert to radians
	const auto lat1 = lat1_p * PI / 180.0;
	const auto lon1 = lon1_p * PI / 180.0;
	const auto lat2 = lat2_p * PI / 180.0;
	const auto lon2 = lon2_p * PI / 180.0;

	const auto dlat = lat2 - lat1;
	const auto dlon = lon2 - lon1;

	const auto a =
	    std::pow(std::sin(dlat / 2.0), 2.0) + std::cos(lat1) * std::cos(lat2) * std::pow(std::sin(dlon / 2.0), 2.0);
	const auto c = 2.0 * std::atan2(std::sqrt(a), std::sqrt(1.0 - a));

	return R * c;
}

} // namespace util

} // namespace sgl
