#include "spatial/util/math.hpp"
#include "duckdb/common/string_util.hpp"

namespace duckdb {

#if SPATIAL_USE_GEOS
// We've got this exposed upstream, we just need to wait for the next release
extern "C" int geos_d2sfixed_buffered_n(double f, uint32_t precision, char *result);

void MathUtil::format_coord(double x, double y, vector<char> &buffer, int32_t precision) {
	D_ASSERT(precision >= 0 && precision <= 15);

	char buf[51];
	auto res_x = geos_d2sfixed_buffered_n(x, 15, buf);
	buf[res_x++] = ' ';
	auto res_y = geos_d2sfixed_buffered_n(y, 15, buf + res_x);
	buffer.insert(buffer.end(), buf, buf + res_x + res_y);
}

void MathUtil::format_coord(double d, vector<char> &buffer, int32_t precision) {
	D_ASSERT(precision >= 0 && precision <= 15);
	char buf[25];
	auto len = geos_d2sfixed_buffered_n(d, 15, buf);
	buffer.insert(buffer.end(), buf, buf + len);
}

string MathUtil::format_coord(double d) {
	char buf[25];
	auto len = geos_d2sfixed_buffered_n(d, 15, buf);
	buf[len] = '\0';
	return string {buf};
}

string MathUtil::format_coord(double x, double y) {
	char buf[51];
	auto res_x = geos_d2sfixed_buffered_n(x, 15, buf);
	buf[res_x++] = ' ';
	auto res_y = geos_d2sfixed_buffered_n(y, 15, buf + res_x);
	buf[res_x + res_y] = '\0';
	return string {buf};
}

string MathUtil::format_coord(double x, double y, double zm) {
	char buf[76];
	auto res_x = geos_d2sfixed_buffered_n(x, 15, buf);
	buf[res_x++] = ' ';
	auto res_y = geos_d2sfixed_buffered_n(y, 15, buf + res_x);
	buf[res_x + res_y++] = ' ';
	auto res_zm = geos_d2sfixed_buffered_n(zm, 15, buf + res_x + res_y);
	buf[res_x + res_y + res_zm] = '\0';
	return string {buf};
}

string MathUtil::format_coord(double x, double y, double z, double m) {
	char buf[101];
	auto res_x = geos_d2sfixed_buffered_n(x, 15, buf);
	buf[res_x++] = ' ';
	auto res_y = geos_d2sfixed_buffered_n(y, 15, buf + res_x);
	buf[res_x + res_y++] = ' ';
	auto res_z = geos_d2sfixed_buffered_n(z, 15, buf + res_x + res_y);
	buf[res_x + res_y + res_z++] = ' ';
	auto res_m = geos_d2sfixed_buffered_n(m, 15, buf + res_x + res_y + res_z);
	buf[res_x + res_y + res_z + res_m] = '\0';
	return string {buf};
}

#else

void MathUtil::format_coord(double x, double y, vector<char> &buffer, int32_t precision) {
	D_ASSERT(precision >= 0 && precision <= 15);
	auto fmt_str = StringUtil::Format("%%.%df %%.%df", precision, precision);
	auto str = StringUtil::Format(fmt_str, x, y);
	buffer.insert(buffer.end(), str.c_str(), str.c_str() + str.size());
}

void MathUtil::format_coord(double d, vector<char> &buffer, int32_t precision) {
	D_ASSERT(precision >= 0 && precision <= 15);
	auto fmt_str = StringUtil::Format("%%.%df", precision);
	auto str = StringUtil::Format(fmt_str, d);
	buffer.insert(buffer.end(), str.c_str(), str.c_str() + str.size());
}

string MathUtil::format_coord(double d) {
	return StringUtil::Format("%G", d);
}

string MathUtil::format_coord(double x, double y) {
	return StringUtil::Format("%G %G", x, y);
}

string MathUtil::format_coord(double x, double y, double zm) {
	return StringUtil::Format("%G %G %G", x, y, zm);
}

string MathUtil::format_coord(double x, double y, double z, double m) {
	return StringUtil::Format("%G %G %G %G", x, y, z, m);
}

#endif

} // namespace duckdb