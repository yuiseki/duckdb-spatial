#pragma once

#ifndef SGL_ASSERT
#ifdef NDEBUG
#define SGL_ASSERT(x) ((void)0)
#else
#include <cassert>
#define SGL_ASSERT(x) assert(x)
#endif
#endif

#include <cstdint>
#include <cstring>
#include <string>
#include <cmath>

namespace sgl {

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
	explicit geometry(const geometry_type type) : type(type) {}

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

	void append_part(geometry *part);

	const uint8_t *get_vertex_data() const;
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
	if(value) {
		flag |= 0x01;
	} else {
		flag &= ~0x01;
	}
	return value;
}

inline bool geometry::set_m(const bool value) {
	if(value) {
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
	if(size == 0) {
		SGL_ASSERT(data == nullptr);
		return nullptr;
	}

	auto part = get_first_part();
	SGL_ASSERT(part != nullptr);

	for(uint32_t i = 0; i < n; i++) {
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
	switch(type) {
		case geometry_type::POINT: return "POINT";
		case geometry_type::LINESTRING: return "LINESTRING";
		case geometry_type::POLYGON: return "POLYGON";
		case geometry_type::MULTI_POINT: return "MULTIPOINT";
		case geometry_type::MULTI_LINESTRING: return "MULTILINESTRING";
		case geometry_type::MULTI_POLYGON: return "MULTIPOLYGON";
		case geometry_type::MULTI_GEOMETRY: return "GEOMETRYCOLLECTION";
		default: return "INVALID";
	}
}

} // namespace sgl

//--------------------------------------------------------------------------
// Operations
//--------------------------------------------------------------------------

namespace sgl {

namespace linestring {
inline geometry make_empty() {
	return geometry(geometry_type::LINESTRING);
}

inline bool is_closed(const geometry *geom) {
	SGL_ASSERT(geom->get_type() == geometry_type::LINESTRING);

	if(geom->get_count() < 2) {
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

	if(count < 3) {
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

	for(uint32_t i = 1; i < count - 1; i++) {
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
	if(count < 2) {
		return 0.0;
	}

	const auto vertex_data = geom->get_vertex_data();
	const auto vertex_size = geom->get_vertex_size();

	auto length = 0.0;

	vertex_xy prev = {0};
	vertex_xy next = {0};

	memcpy(&prev, vertex_data, sizeof(vertex_xy));

	for(uint32_t i = 1; i < count; i++) {
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

inline geometry make_empty() {
	return geometry(geometry_type::POLYGON);
}

inline double area(const geometry *geom) {
	SGL_ASSERT(geom->get_type() == geometry_type::POLYGON);
	double area = 0.0;

	auto part = geom->get_first_part();
	if(!part) {
		return area;
	}

	area += std::abs(linestring::signed_area(part));

	while(part != geom->get_last_part()) {
		part = part->get_next();
		area -= std::abs(linestring::signed_area(part));
	}

	return area;
}

} // namespace polygon

namespace multi_point {

inline geometry make_empty() {
	return geometry(geometry_type::MULTI_POINT);
}

};

namespace multi_line_string {

inline geometry make_empty() {
	return geometry(geometry_type::MULTI_LINESTRING);
}

};

namespace multi_polygon {
inline geometry make_empty() {
	return geometry(geometry_type::MULTI_POLYGON);
}

inline double area(const geometry *geom) {
	SGL_ASSERT(geom->get_type() == geometry_type::MULTI_POLYGON);

	const auto tail = geom->get_last_part();
	if(!tail) {
		return 0.0;
	}

	double area = 0.0;
	auto part = tail;
	do {
		part = part->get_next();
		area += polygon::area(part);
	} while(part != tail);

	return area;
}

}

namespace multi_geometry {

inline geometry make_empty() {
	return geometry(geometry_type::MULTI_GEOMETRY);
}

inline double area(const geometry *geom);

};

namespace ops {

double area(const geometry *geom);

} // namespace ops

} // namespace sgl

//--------------------------------------------------------------------------
// Implementation
//--------------------------------------------------------------------------

namespace sgl {
namespace ops {

inline double area(const geometry *geom) {
	switch(geom->get_type()) {
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

} // namespace ops
} // namespace sgl

namespace sgl {
namespace multi_geometry {

inline double area(const geometry *geom) {
	SGL_ASSERT(geom->get_type() == geometry_type::MULTI_GEOMETRY);

	const auto tail = geom->get_last_part();
	if(!tail) {
		return 0.0;
	}

	double area = 0.0;
	auto part = tail;
	do {
		part = part->get_next();
		area += ops::area(part);
	} while(part != tail);

	return area;
}

}
}
