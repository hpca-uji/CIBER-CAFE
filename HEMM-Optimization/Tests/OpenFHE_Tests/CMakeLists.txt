cmake_minimum_required(VERSION 3.15)

# Option to force using g++/gcc instead of the detected compiler (useful if clang causes issues)
option(FORCE_GXX "Force selection of g++/gcc as the C/C++ compiler" OFF)
if(FORCE_GXX)
  find_program(GPP_PATH g++)
  find_program(GCC_PATH gcc)
  if(GPP_PATH)
    message(STATUS "FORCE_GXX=ON: forcing C++ compiler to ${GPP_PATH}")
    set(CMAKE_CXX_COMPILER "${GPP_PATH}" CACHE PATH "C++ compiler" FORCE)
  else()
    message(WARNING "FORCE_GXX=ON but g++ not found in PATH. Ignoring.")
  endif()
  if(GCC_PATH)
    message(STATUS "FORCE_GXX=ON: forcing C compiler to ${GCC_PATH}")
    set(CMAKE_C_COMPILER "${GCC_PATH}" CACHE PATH "C compiler" FORCE)
  endif()
endif()

project(matMultTest C CXX)

set(CMAKE_CXX_STANDARD 17)
set(CXX_STANDARD_REQUIRED ON)

# Opciones de compilación para controlar debug y uso de pragmas OpenMP
option(ENABLE_DEBUG "Enable debug prints in matMultTest" OFF)
option(USE_OMP_TASKLOOP "Enable OpenMP taskloop/taskgroup pragmas (may be unstable on some compilers)" OFF)

# Detectar clang para desactivar taskloop por defecto si la versión es antigua
if(CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
  if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS 15.0)
    if(NOT DEFINED USE_OMP_TASKLOOP)
      message(WARNING "Clang < 15 detected: defaulting USE_OMP_TASKLOOP=OFF to avoid frontend crashes.")
      set(USE_OMP_TASKLOOP OFF CACHE BOOL "Disable OpenMP taskloop for clang < 15" FORCE)
    endif()
  endif()
endif()

# If the user did not explicitly request USE_OMP_TASKLOOP, try to detect compiler support
if(NOT DEFINED USE_OMP_TASKLOOP)
  message(STATUS "Detecting OpenMP taskloop support...")
  include(CheckCXXSourceCompiles)
  set(CMAKE_REQUIRED_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  check_cxx_source_compiles(
    "#include <omp.h>\nint main(){\n#pragma omp parallel\n{\n#pragma omp taskloop for\n for(int i=0;i<1;i++);\n}\n return 0;\n}"
    HAVE_OPENMP_TASKLOOP)
  if(HAVE_OPENMP_TASKLOOP)
    message(STATUS "Compiler supports OpenMP taskloop/taskgroup. Setting USE_OMP_TASKLOOP=ON by default.")
    set(USE_OMP_TASKLOOP ON CACHE BOOL "Enable OpenMP taskloop if supported" FORCE)
  else()
    message(STATUS "Compiler does NOT support OpenMP taskloop/taskgroup. Setting USE_OMP_TASKLOOP=OFF by default.")
    set(USE_OMP_TASKLOOP OFF CACHE BOOL "Disable OpenMP taskloop if unsupported" FORCE)
  endif()
endif()


# Expect user to pass -DOPENFHE_ROOT=/path/to/openfhe_install
if(NOT DEFINED OPENFHE_ROOT)
    message(FATAL_ERROR "You must provide -DOPENFHE_ROOT=/path/to/openfhe_install")
endif()

set(OPENFHE_INCLUDE_DIR ${OPENFHE_ROOT}/include)
set(OPENFHE_LIB_DIR ${OPENFHE_ROOT}/lib)

# Add includes that match the installation layout
include_directories(
    ${OPENFHE_INCLUDE_DIR}
    ${OPENFHE_INCLUDE_DIR}/openfhe
    ${OPENFHE_INCLUDE_DIR}/openfhe/pke
    ${OPENFHE_INCLUDE_DIR}/openfhe/core
    ${OPENFHE_INCLUDE_DIR}/openfhe/binfhe
)

## Añadir utils.cpp desde la carpeta padre (../utils.cpp) como fuente compartida
add_executable(matMultTest matMultTest.cpp)

## Asegurar que el compilador busque headers en el directorio fuente del proyecto y en la carpeta padre
target_include_directories(matMultTest PRIVATE ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/..)

# Convertir opciones en macros para el compilador (mover después de add_executable)
if(ENABLE_DEBUG)
  target_compile_definitions(matMultTest PRIVATE ENABLE_DEBUG=1)
else()
  target_compile_definitions(matMultTest PRIVATE ENABLE_DEBUG=0)
endif()

if(USE_OMP_TASKLOOP)
  target_compile_definitions(matMultTest PRIVATE USE_OMP_TASKLOOP=1)
else()
  target_compile_definitions(matMultTest PRIVATE USE_OMP_TASKLOOP=0)
endif()

# Report chosen options to the user
message(STATUS "CMake option: ENABLE_DEBUG=${ENABLE_DEBUG}")
message(STATUS "CMake option: USE_OMP_TASKLOOP=${USE_OMP_TASKLOOP}")

# Extra warning if taskloop is enabled but compiler is known to be problematic
if(CMAKE_CXX_COMPILER_ID STREQUAL "Clang" AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS 15.0 AND USE_OMP_TASKLOOP)
  message(WARNING "USE_OMP_TASKLOOP=ON but Clang version < 15 detected: OpenMP taskloop/taskgroup may trigger frontend crashes. Consider disabling USE_OMP_TASKLOOP or using a newer compiler.")
endif()

# Link the specific libraries present in this OpenFHE layout
find_library(OPENFHE_PKE NAMES OPENFHEpke PATHS ${OPENFHE_LIB_DIR})
find_library(OPENFHE_CORE NAMES OPENFHEcore PATHS ${OPENFHE_LIB_DIR})
find_library(OPENFHE_BINFHE NAMES OPENFHEbinfhe PATHS ${OPENFHE_LIB_DIR})

# Find dependencies required by OpenFHE installation (NTL, GMP)
find_library(NTL NAMES ntl NAMES ntl.so ntl PATHS ${OPENFHE_LIB_DIR})
find_library(GMP NAMES gmp PATHS ${OPENFHE_LIB_DIR})
find_library(HEXl NAMES hexl hexl.so libhexl.so hexl PATHS ${OPENFHE_LIB_DIR})

if(NOT OPENFHE_PKE OR NOT OPENFHE_CORE)
    message(FATAL_ERROR "Could not find required OPENFHE libraries in ${OPENFHE_LIB_DIR}")
endif()

target_link_libraries(matMultTest PRIVATE ${OPENFHE_PKE} ${OPENFHE_CORE} ${OPENFHE_BINFHE} ${NTL} ${GMP} ${HEXl})

# Option to enable/disable OpenMP at configure time
option(ENABLE_OPENMP "Enable use of OpenMP (linking and flags)" ON)

if(ENABLE_OPENMP)
  # Try to find an OpenMP installation; do not FATAL if missing so user can disable via flag
  find_package(OpenMP)
  if(OpenMP_CXX_FOUND)
    message(STATUS "OpenMP found: enabling OpenMP linking and flags")
    target_link_libraries(matMultTest PRIVATE OpenMP::OpenMP_CXX)
    target_compile_definitions(matMultTest PRIVATE ENABLE_OPENMP=1)
  else()
    message(WARNING "ENABLE_OPENMP=ON but OpenMP not found. Building without OpenMP.")
    target_compile_definitions(matMultTest PRIVATE ENABLE_OPENMP=0)
  endif()
else()
  message(STATUS "OpenMP disabled via ENABLE_OPENMP=OFF")
  target_compile_definitions(matMultTest PRIVATE ENABLE_OPENMP=0)
endif()

# Set rpath so executable finds the libs at runtime
set_target_properties(matMultTest PROPERTIES
    BUILD_RPATH "${OPENFHE_LIB_DIR}"
    INSTALL_RPATH "${OPENFHE_LIB_DIR}"
)

# C++ standard
target_compile_features(matMultTest PRIVATE cxx_std_17)

message(STATUS "OPENFHE_INCLUDE_DIR=${OPENFHE_INCLUDE_DIR}")
message(STATUS "OPENFHE_LIB_DIR=${OPENFHE_LIB_DIR}")
