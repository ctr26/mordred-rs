cmake_minimum_required(VERSION 3.15)
project(mordred-ffi LANGUAGES C CXX)

# --------------------------------------------------------------------------
# Configuration
# --------------------------------------------------------------------------
set(MORDRED_RS_ROOT "${CMAKE_CURRENT_SOURCE_DIR}/../.." CACHE PATH
    "Root of the mordred-rs workspace")
set(MORDRED_BUILD_TYPE "release" CACHE STRING
    "Cargo build profile: 'debug' or 'release'")

if(MORDRED_BUILD_TYPE STREQUAL "release")
  set(CARGO_PROFILE_FLAG "--release")
  set(CARGO_TARGET_SUBDIR "release")
else()
  set(CARGO_PROFILE_FLAG "")
  set(CARGO_TARGET_SUBDIR "debug")
endif()

set(MORDRED_TARGET_DIR "${MORDRED_RS_ROOT}/target/${CARGO_TARGET_SUBDIR}")
set(MORDRED_INCLUDE_DIR "${MORDRED_RS_ROOT}/mordred-ffi/include")

# --------------------------------------------------------------------------
# Build the Rust crate via cargo
# --------------------------------------------------------------------------
add_custom_target(mordred_cargo ALL
  COMMAND cargo build ${CARGO_PROFILE_FLAG} -p mordred-ffi
  WORKING_DIRECTORY "${MORDRED_RS_ROOT}"
  COMMENT "Building mordred-ffi Rust crate"
)

# --------------------------------------------------------------------------
# Imported library target: mordred::mordred
# --------------------------------------------------------------------------
# Prefer the static library; fall back to the shared library.
if(WIN32)
  set(_STATIC_LIB "${MORDRED_TARGET_DIR}/mordred_ffi.lib")
  set(_SHARED_LIB "${MORDRED_TARGET_DIR}/mordred_ffi.dll")
else()
  set(_STATIC_LIB "${MORDRED_TARGET_DIR}/libmordred_ffi.a")
  set(_SHARED_LIB "${MORDRED_TARGET_DIR}/libmordred_ffi.so")
  if(APPLE)
    set(_SHARED_LIB "${MORDRED_TARGET_DIR}/libmordred_ffi.dylib")
  endif()
endif()

if(EXISTS "${_STATIC_LIB}" OR NOT EXISTS "${_SHARED_LIB}")
  add_library(mordred_mordred STATIC IMPORTED GLOBAL)
  set_target_properties(mordred_mordred PROPERTIES
    IMPORTED_LOCATION "${_STATIC_LIB}"
  )
else()
  add_library(mordred_mordred SHARED IMPORTED GLOBAL)
  set_target_properties(mordred_mordred PROPERTIES
    IMPORTED_LOCATION "${_SHARED_LIB}"
  )
endif()

add_dependencies(mordred_mordred mordred_cargo)

target_include_directories(mordred_mordred INTERFACE "${MORDRED_INCLUDE_DIR}")

# Platform-specific link dependencies for the static library.
if(UNIX AND NOT APPLE)
  set_property(TARGET mordred_mordred APPEND PROPERTY
    INTERFACE_LINK_LIBRARIES dl pthread m)
elseif(APPLE)
  set_property(TARGET mordred_mordred APPEND PROPERTY
    INTERFACE_LINK_LIBRARIES dl pthread "-framework Security")
elseif(WIN32)
  set_property(TARGET mordred_mordred APPEND PROPERTY
    INTERFACE_LINK_LIBRARIES ws2_32 userenv bcrypt ntdll)
endif()

# Provide a namespaced alias.
add_library(mordred::mordred ALIAS mordred_mordred)
