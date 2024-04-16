#ifndef LICOM3_KOKKOS_SRC_HEAD_CPP_POP_GRID_HORZ_MOD_H_
#define LICOM3_KOKKOS_SRC_HEAD_CPP_POP_GRID_HORZ_MOD_H_

namespace CppPOPGridHorzMod {

constexpr const char POP_FIELD_KIND_UNKNOWN[]            = "unknown";
constexpr const char POP_FIELD_KIND_SCALAR[]             = "scalar";
constexpr const char POP_FIELD_KIND_VECTOR[]             = "vector";
constexpr const char POP_FIELD_KIND_ANGLE[]              = "angle";
constexpr const char POP_FIELD_KIND_NO_UPDATE[]          = "noUpdate";

constexpr const char POP_GRID_HORZ_KIND_UNKNOWN[]        = "unknow";
constexpr const char POP_GRID_HORZ_KIND_LAT_LON[]        = "LatLon";
constexpr const char POP_GRID_HORZ_KIND_DISPLACED_POLE[] = "DisplacedPole";
constexpr const char POP_GRID_HORZ_KIND_TRIPOLE[]        = "Tripole";
constexpr const char POP_GRID_HORZ_KIND_GEODESIC[]       = "Geodesic";

constexpr const char POP_GRID_HORZ_STAGGER_UNKNOW[]      = "Unknow";
constexpr const char POP_GRID_HORZ_STAGGER_A[]           = "A";
constexpr const char POP_GRID_HORZ_STAGGER_B_NE[]        = "B_N";
constexpr const char POP_GRID_HORZ_STAGGER_C_NE[]        = "C_N";
constexpr const char POP_GRID_HORZ_STAGGER_Z[]           = "Z";

constexpr const char POP_GRID_HORZ_LOC_UNKNOW[]          = "Unknow";
constexpr const char POP_GRID_HORZ_LOC_CENTER[]          = "Center";
constexpr const char POP_GRID_HORZ_LOC_N_FACE[]          = "NFace";
constexpr const char POP_GRID_HORZ_LOC_S_FACE[]          = "SFace";
constexpr const char POP_GRID_HORZ_LOC_E_FACE[]          = "EFace";
constexpr const char POP_GRID_HORZ_LOC_W_FACE[]          = "WFace";
constexpr const char POP_GRID_HORZ_LOC_NE_CORNER[]       = "NECorner";
constexpr const char POP_GRID_HORZ_LOC_NW_CORNER[]       = "NWCorner";
constexpr const char POP_GRID_HORZ_LOC_SE_CORNER[]       = "SECorner";
constexpr const char POP_GRID_HORZ_LOC_SW_CORNER[]       = "SWCorner";
constexpr const char POP_GRID_HORZ_LOC_NO_UPDATE[]       = "NoUpdate";

// "unknown"
constexpr int FLAG_POP_FIELD_KIND_UNKNOWN              = 1;
// "scalar"
constexpr int FLAG_POP_FIELD_KIND_SCALAR               = 2;
// "vector"
constexpr int FLAG_POP_FIELD_KIND_VECTOR               = 3;
// "angle"
constexpr int FLAG_POP_FIELD_KIND_ANGLE                = 4;
// "noUpdate"
constexpr int FLAG_POP_FIELD_KIND_NO_UPDATE            = 5;

// "unknow"
constexpr int FLAG_POP_GRID_HORZ_KIND_UNKNOWN          = 6;
// "LatLon"
constexpr int FLAG_POP_GRID_HORZ_KIND_LAT_LON          = 7;
// "DisplacedPole"
constexpr int FLAG_POP_GRID_HORZ_KIND_DISPLACED_POLE   = 8;
// "Tripole"
constexpr int FLAG_POP_GRID_HORZ_KIND_TRIPOLE          = 9;
// "Geodesic"
constexpr int FLAG_POP_GRID_HORZ_KIND_GEODESIC         = 10;

// "Unknow"
constexpr int FLAG_POP_GRID_HORZ_STAGGER_UNKNOW        = 11;
// "A"
constexpr int FLAG_POP_GRID_HORZ_STAGGER_A             = 12;
// "B_N"
constexpr int FLAG_POP_GRID_HORZ_STAGGER_B_NE          = 13;
// "C_N"
constexpr int FLAG_POP_GRID_HORZ_STAGGER_C_NE          = 14;
// "Z"
constexpr int FLAG_POP_GRID_HORZ_STAGGER_Z             = 15;

// "Unknow"
constexpr int FLAG_POP_GRID_HORZ_LOC_UNKNOW            = 16;
// "Center"
constexpr int FLAG_POP_GRID_HORZ_LOC_CENTER            = 17;
// "NFace"
constexpr int FLAG_POP_GRID_HORZ_LOC_N_FACE            = 18;
// "SFace"
constexpr int FLAG_POP_GRID_HORZ_LOC_S_FACE            = 19;
// "EFace"
constexpr int FLAG_POP_GRID_HORZ_LOC_E_FACE            = 20;
// "WFace"
constexpr int FLAG_POP_GRID_HORZ_LOC_W_FACE            = 21;
// "NECorner"
constexpr int FLAG_POP_GRID_HORZ_LOC_NE_CORNER         = 22;
// "NWCorner"
constexpr int FLAG_POP_GRID_HORZ_LOC_NW_CORNER         = 23;
// "SECorner"
constexpr int FLAG_POP_GRID_HORZ_LOC_SE_CORNER         = 24;
// "SWCorner"
constexpr int FLAG_POP_GRID_HORZ_LOC_SW_CORNER         = 25;
// "NoUpdate"
constexpr int FLAG_POP_GRID_HORZ_LOC_NO_UPDATE         = 26;

} // namespace CppPOPGridHorzMod

#endif // LICOM3_KOKKOS_SRC_HEAD_CPP_POP_GRID_HORZ_MOD_H_
