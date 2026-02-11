#!/bin/bash
#SBATCH --job-name=HEMM-Optimization
#SBATCH --exclusive
#SBATCH --time=7-00:00:00
#SBATCH --output=HEMM-Optimization_%j.out

# Configuración
libInstall=~/.local
raizProyecto=$(pwd)

# Configuración específica por nodo/partición
if [ -n "$SLURM_JOB_PARTITION" ] && [ -n "$SLURMD_NODENAME" ]; then
    # Directorio específico para cada nodo en el clúster
    libInstall=$libInstall/${SLURM_JOB_PARTITION}.${SLURMD_NODENAME}
elif [ -n "$HOSTNAME" ]; then
    # Fallback: usar hostname si no hay SLURM
    libInstall=$libInstall/$(hostname -s)
fi

# Crear directorio de instalación si no existe
mkdir -p "$libInstall"

gmpVers="6.2.1"
gmpRuta=~/gmp-$gmpVers

gf2x=False                      # Cambiar a True para instalar GF2X
gf2xVers="1.3.0"
gf2xRuta=~/gf2x-$gf2xVers

ntlVers="11.5.1"
ntlRuta=~/ntl-$ntlVers

hexlBranch="main"
hexlRuta=~/.git/HEXL

openFHEBranch=""
openFHERuta=~/.git/OpenFHEConfig
openFHEArgs=("-DNATIVE_SIZE=64" "-DNATIVE_SIZE=128")

helibBranch=""
helibRuta=~/.git/HElib
helibArgs=("")

sealRuta=~/.git/SEAL
sealArgs=("")

formato_fecha() {
    date '+%Y%m%d_%H%M%S%2N'
}

ruta_job=$raizProyecto/Resultados/HEMM-Optimization_
if [ -n "$SLURM_JOB_PARTITION" ] && [ -n "$SLURMD_NODENAME" ]; then
    # Directorio específico para cada nodo en el clúster
    ruta_job=$ruta_job${SLURM_JOB_PARTITION}.${SLURMD_NODENAME}_
elif [ -n "$HOSTNAME" ]; then
    # Fallback: usar hostname si no hay SLURM
    ruta_job=$ruta_job$(hostname -s)_
fi
ruta_job=$ruta_job$(formato_fecha)
logCombinado=$ruta_job/ejecucion.log
ficResul=results.csv

# Opciones
copiarfichResult=False  # Cambiar a True para copiar el fichero de resultados si existe
colorOutput=False       # Habilitar salida con colores
check=True              # Cambiar a True para hacer 'make check'
remove=True             # Cambiar a True para eliminar los archivos fuente
numNegativos=False      # Permitir números negativos en las matrices generadas
ficMat1=mat1            # Nombre del archivo de matriz 1
ficMat2=mat2            # Nombre del archivo de matriz 2
precision=16            # Precisión de los números en las matrices (número de dígitos en la parte fraccionaria)
timeoutSec=3600         # 60 minutos

#Detalles de los tests
openmpInAlg=OFF  # Habilitar OpenMP en los Algoritmos
repeticiones=10         # Número de repeticiones de cada test
reiteraciones=(1)    # Número de reiteraciones de la multiplicación
tamsMatriz=(2 3 4 6 8 12 16 24 32 48 64 96 128 192 256 384 512 768 1024 1536 2048 3072 4096) # Tamaños de las matrices a probar
algRectangulares=("Traditional" "RowsXCols" "HaleviShoup")  # Algoritmos para matrices rectangulares

OpenFHEAlg=("HEMatMult" "RizomiliotisTriakosia" "RowsXCols" "HaleviShoup" "Traditional" "Strassen1x1x1" "StrassenRizomiliotisTriakosia")
OpenFHERingDim=(8192 16384 32768 65536 131072 262144 524288 1048576 2097152 4194304 8388608 16777216 33554432)
OpenFHEMultDepth=(1 2 3 4 5 6 7 8 9 10)

HElibAlg=("HEMatMult" "RizomiliotisTriakosia" "RowsXCols" "HaleviShoup" "Traditional" "Strassen1x1x1" "StrassenRizomiliotisTriakosia")
HElibParams=("16384 119 2" "32768 358 6" "32768 299 3" "32768 239 2" "65536 725 8" "65536 717 6" "65536 669 4" "65536 613 3" "65536 558 2" "131072 1445 8" "131072 1435 6" "131072 1387 5" "131072 1329 4" "131072 1255 3" "131072 1098 2" "262144 2940 8" "262144 2870 6" "262144 2763 5" "262144 2646 4" "262144 2511 3" "262144 2234 2")
HElibPrecisions=(40)

# Crear directorio de trabajo primero
mkdir -p $ruta_job

if [ "$colorOutput" = "True" ]; then
    RED="\e[31m"
    GREEN="\e[32m"
    YELLOW="\e[33m"
    ORANGE="\e[33;1m"
    BLUE="\e[34m"
    BLACK="\e[30m"
    RESET="\e[0m"
else
    RED=""
    GREEN=""
    YELLOW=""
    ORANGE=""
    BLUE=""
    BLACK=""
    RESET=""
fi

# Funciones para logging con timestamp y prefijos
log_info() {
    echo "${BLUE}[$(formato_fecha)] ℹ️ [INFO]    $1${RESET}" | tee -a $logCombinado
}
log_error() {
    echo "${RED}[$(formato_fecha)] ❌ [ERROR]   $1${RESET}" | tee -a $logCombinado >&2
}
log_command() {
    echo "${BLACK}[$(formato_fecha)] 💻 [CMD]     $1${RESET}" | tee -a $logCombinado
}
log_test() {
    echo "${ORANGE}[$(formato_fecha)] 🧪 [TEST]    $1${RESET}" | tee -a $logCombinado
}
log_install() {
    echo "${YELLOW}[$(formato_fecha)] 🛠️ [INSTALL] $1${RESET}" | tee -a $logCombinado
}
log_result() {
    echo "${GREEN}[$(formato_fecha)] 📊 [RESULT]  $1${RESET}" | tee -a $logCombinado
}

# Función para capturar salida de comandos con formato
log_output() {
    while IFS= read -r line; do
        echo "${GREEN}[$(formato_fecha)] 📤 [OUT]     $line${RESET}" | tee -a $logCombinado
    done
}
log_install_output() {
    while IFS= read -r line; do
        echo "${YELLOW}[$(formato_fecha)] 🛠️ [INSTALL] $line${RESET}" | tee -a $logCombinado
    done
}

log_error_output() {
    while IFS= read -r line; do
        echo "${RED}[$(formato_fecha)] ❌ [ERROR]   $line${RESET}" | tee -a $logCombinado
    done
}
log_test_output() {
    while IFS= read -r line; do
        echo "${ORANGE}[$(formato_fecha)] 🧪 [TEST]    $line${RESET}" | tee -a $logCombinado
    done
}

# Función para mostrar información del sistema de forma elegante
mostrar_info_sistema() {
    log_info "🔍 Detectando información del sistema..."
    log_info "┌─────────────────────────────────────────────────────────────┐"
    log_info "│                    INFORMACIÓN DEL SISTEMA                  │"
    log_info "├─────────────────────────────────────────────────────────────┤"
    
    log_info "│ 📅 Fecha: $(date '+%Y-%m-%d %H:%M:%S')"
    log_info "│ 👤 Usuario: $(whoami)"
    log_info "│ 🖥️  Sistema: $(uname -s) $(uname -r)"
    if [ -n "$SLURM_JOB_PARTITION" ] && [ -n "$SLURMD_NODENAME" ]; then
        log_info "│ 🖥️  Nodo SLURM: $SLURMD_NODENAME (Partición: $SLURM_JOB_PARTITION)"
    else
        log_info "│ 🖥️  Nodo: $(hostname -s)"
    fi
    
    # CPU Information
    log_info "│ 🖥️  PROCESADOR:"
    local cpu_model=$(grep "model name" /proc/cpuinfo | head -1 | cut -d: -f2 | sed 's/^[ \t]*//')
    local cpu_cores=$(grep "cpu cores" /proc/cpuinfo | head -1 | cut -d: -f2 | sed 's/^[ \t]*//')
    local physical_cpus=$(grep "physical id" /proc/cpuinfo | sort | uniq | wc -l)
    log_info "│   Modelo: $cpu_model"
    log_info "│   CPUs físicas: $physical_cpus | Cores por CPU: $cpu_cores"
    log_info "│   Threads totales: $(nproc)"
    
    # Memory Information  
    log_info "│ 💾 MEMORIA:"
    local mem_info=$(free -h | grep "Mem:")
    local total_mem=$(echo $mem_info | awk '{print $2}')
    local available_mem=$(echo $mem_info | awk '{print $7}')
    log_info "│   Total: $total_mem | Disponible: $available_mem"
    
    # Operating System
    log_info "│ 🐧 SISTEMA OPERATIVO:"
    if [ -f /etc/os-release ]; then
        local os_name=$(grep "PRETTY_NAME" /etc/os-release | cut -d'"' -f2)
        log_info "│   $os_name"
        log_info "│   Kernel: $(uname -r)"
    fi
    
    # Development Tools
    log_info "│ 🛠️  HERRAMIENTAS DE DESARROLLO:"
    # CMake
    if command -v cmake >/dev/null 2>&1; then
        local cmake_ver=$(cmake --version | head -1 | awk '{print $3}')
        log_info "│   CMake: v$cmake_ver"
    else
        log_info "│   CMake: ❌ No instalado"
    fi
    # GCC
    if command -v gcc >/dev/null 2>&1; then
        local gcc_ver=$(gcc --version | head -1 | awk '{print $4}')
        log_info "│   GCC: v$gcc_ver"
    else
        log_info "│   GCC: ❌ No instalado"
    fi
    # G++
    if command -v g++ >/dev/null 2>&1; then
        local gpp_ver=$(g++ --version | head -1 | awk '{print $4}')
        log_info "│   G++: v$gpp_ver"
    else
        log_info "│   G++: ❌ No instalado"
    fi
    # Clang
    if command -v clang++ >/dev/null 2>&1; then
        local clang_ver=$(clang++ --version | head -1 | awk '{print $4}')
        log_info "│   Clang++: v$clang_ver"
    else
        log_info "│   Clang++: ❌ No instalado"
    fi
    
    log_info "├────────────────────────────────────────────────────────────────────┤"
    log_info "│  🚀 Iniciando instalación de dependencias para HEMM-Optimization   │"
    log_info "├────────────────────────────────────────────────────────────────────┤"
    log_info "│ 📁 Directorio de instalación: $libInstall"
    log_info "│ 📋 Configuración:"
    log_info "│   • GMP v$gmpVers"
    log_info "│   • GF2X v$gf2xVers (habilitado: $gf2x)"
    log_info "│   • NTL v$ntlVers"
    log_info "│   • OpenFHE (configuraciones: ${openFHEArgs[*]})"
    log_info "│   • Ejecutar tests: $check"
    log_info "│   • Limpiar archivos: $remove"
    log_info "│ 📝 Log completo: $logCombinado"
    log_info "└────────────────────────────────────────────────────────────────────┘"
}

# Información inicial del proceso
mostrar_info_sistema

# Copiar o crear ficheros de matrices
maxTam=0
for tam in "${tamsMatriz[@]}"; do
    if (( tam > maxTam )); then
        maxTam=$tam
    fi
done
ficMat1="$raizProyecto/${ficMat1}_${maxTam}.csv"
ficMat2="$raizProyecto/${ficMat2}_${maxTam}.csv"

# Helper: genera una cadena numérica aleatoria de 'precision' dígitos
# Generador rápido de matrices CSV usando awk
# - Usa /dev/urandom y awk para generar líneas con números 0.<dígitos>
# - FALLBACK_LEGACY_GENERATOR=1 fuerza el uso del generador bash original
generar_matriz_csv() {
    local tam=$1
    # Método rápido con awk + /dev/urandom
    # awk no tiene acceso directo a /dev/urandom en todas las plataformas, así que
    # creamos una función generadora en awk que usa srand() y rand() y formatea
    # la parte fraccionaria con la longitud indicada.
    # Pasamos numNegativos a awk como la variable NEG para evaluar correctamente
    awk -v N="$tam" -v P="$precision" -v NEG="$numNegativos" 'BEGIN {
        srand();
        for (i=0;i<N;i++){
            line = "";
            for (j=0;j<N;j++){
                # Construir la parte fraccionaria concatenando P dígitos [0-9]
                s = "";
                for (k=0;k<P;k++) {
                    d = int(rand()*10);
                    s = s sprintf("%d", d);
                }
                token = "0." s;
                if (NEG == "True" && rand() < 0.5) {
                    token = "-" token;
                }
                if (j==0) line = token; else line = line "," token;
            }
            print line;
        }
    }'
}

if [ ! -f "$ficMat1" ]; then
    log_info "Fichero de matriz 1 no encontrado, creando nuevo fichero con tamaño máximo $maxTam x $maxTam"
    generar_matriz_csv $maxTam > $ficMat1
    log_info "Fichero $ficMat1 creado correctamente"
fi
if [ ! -f "$ficMat2" ]; then
    log_info "Fichero de matriz 2 no encontrado, creando nuevo fichero con tamaño máximo $maxTam x $maxTam"
    generar_matriz_csv $maxTam > $ficMat2
    log_info "Fichero $ficMat2 creado correctamente"
fi

if [ "$copiarfichResult" = "True" ]; then
    # Activar nullglob: si no hay coincidencias, el glob desaparece (no queda sin expandir)
    shopt -s nullglob
    files=( "$raizProyecto/Resultados/"*"$ficResul" )
    if (( ${#files[@]} > 0 )); then
        log_info "Copiando ${#files[@]} fichero(s) que terminan en $ficResul al directorio de trabajo"
        log_command cp -t "$ruta_job/" -- "${files[@]}"
        cp -t "$ruta_job/" -- "${files[@]}"
        log_info "Fichero(s) copiado(s) correctamente"
    else
        log_info "Fichero de resultados $ficResul no encontrado en el proyecto, se creará uno nuevo en el directorio de trabajo"
    fi
    # Opcional: desactivar para no afectar otras partes del script
    shopt -u nullglob
fi

# Instalación de GMP
if [ ! -f "$libInstall/lib/libgmp.a" ]; then
    log_install "Iniciando instalación de GMP v$gmpVers"
    mkdir -p $gmpRuta || exit
    log_command "cd $gmpRuta"
    cd $gmpRuta || exit
    if [ ! -f "gmp-$gmpVers.tar.xz" ]; then
        log_command "wget https://ftp.gnu.org/gnu/gmp/gmp-$gmpVers.tar.xz"
        wget https://ftp.gnu.org/gnu/gmp/gmp-$gmpVers.tar.xz 2>&1 | log_install_output || exit
    fi
    log_command "tar -xf gmp-$gmpVers.tar.xz"
    tar -xf gmp-$gmpVers.tar.xz 2>&1 | log_install_output || exit
    ( \
        log_command "cd ./gmp-$gmpVers || exit;"; \
        cd ./gmp-$gmpVers || exit; \
        log_command "./configure --prefix=$libInstall"; \
        ./configure --prefix=$libInstall 2>&1 | log_install_output; \
        log_command "make"; \
        make 2>&1 | log_install_output; \
        if [ "$check" = "True" ]; then
            log_test "make check"; \
            make check 2>&1 | log_install_output; \
        fi
        log_install "make install"; \
        make install 2>&1 | log_install_output; \
        log_command "cd .."; \
        cd ..; \
    )
    log_command "rm -rf gmp-$gmpVers"
    rm -rf gmp-$gmpVers || exit
    if [ "$remove" = "True" ]; then
        log_command "rm -rf $gmpRuta"
        rm -rf $gmpRuta || exit
    fi
    log_install "✅ GMP v$gmpVers instalado correctamente"
else
    log_install "⏭️  GMP ya está instalado, saltando..."
fi

# Instalación de GF2X (opcional)
if [ "$gf2x" = "True" ]; then
    if [ ! -f "$libInstall/lib/libgf2x.a" ]; then
        log_install "Iniciando instalación de GF2X v$gf2xVers"
        mkdir -p $gf2xRuta || exit
        log_command "cd $gf2xRuta"
        cd $gf2xRuta || exit
        if [ ! -f "gf2x-$gf2xVers.tar.gz" ]; then
            log_command "wget https://gitlab.inria.fr/-/project/17662/uploads/c46b1047ba841c20d1225ae73ad6e4cd/gf2x-$gf2xVers.tar.gz"
            wget https://gitlab.inria.fr/-/project/17662/uploads/c46b1047ba841c20d1225ae73ad6e4cd/gf2x-$gf2xVers.tar.gz 2>&1 | log_install_output || exit
        fi
        log_command "tar -xf gf2x-$gf2xVers.tar.gz"
        tar -xf gf2x-$gf2xVers.tar.gz 2>&1 | log_install_output || exit
        ( \
            log_command "cd ./gf2x-$gf2xVers || exit;"; \
            cd ./gf2x-$gf2xVers || exit; \
            log_command "./configure --prefix=$libInstall"; \
            ./configure --prefix=$libInstall 2>&1 | log_install_output; \
            log_command "make"; \
            make 2>&1 | log_install_output; \
            if [ "$check" = "True" ]; then
                log_test "make check"; \
                make check 2>&1 | log_install_output; \
            fi
            log_install "make install"; \
            make install 2>&1 | log_install_output; \
            log_command "cd .."; \
            cd ..; \
        )
        log_command "rm -rf gf2x-$gf2xVers"
        rm -rf gf2x-$gf2xVers || exit
        if [ "$remove" = "True" ]; then
            log_command "rm -rf $gf2xRuta"
            rm -rf $gf2xRuta || exit
        fi
        log_install "✅ GF2X v$gf2xVers instalado correctamente"
    else
        log_install "⏭️  GF2X ya está instalado, saltando..."
    fi
    ntlargs="NTL_GF2X_LIB=on GF2X_PREFIX=$libInstall"
else
    log_install "⏭️  GF2X deshabilitado (gf2x=False)"
    ntlargs=""
fi

# Instalación de NTL
if [ ! -f "$libInstall/lib/libntl.a" ]; then
    log_install "Iniciando instalación de NTL v$ntlVers"
    mkdir -p $ntlRuta || exit
    log_command "cd $ntlRuta"
    cd $ntlRuta || exit
    if [ ! -f "ntl-$ntlVers.tar.gz" ]; then
        log_command "wget https://libntl.org/ntl-$ntlVers.tar.gz"
        wget https://libntl.org/ntl-$ntlVers.tar.gz 2>&1 | log_install_output || exit
    fi
    log_command "tar -xf ntl-$ntlVers.tar.gz"
    tar -xf ntl-$ntlVers.tar.gz 2>&1 | log_install_output || exit
    ( \
        log_command "cd ./ntl-$ntlVers/src;"; \
        cd ./ntl-$ntlVers/src; \
        log_command "./configure NTL_EXCEPTIONS=on SHARED=on NTL_STD_CXX11=on NTL_SAFE_VECTORS=off TUNE=generic PREFIX=$libInstall GMP_PREFIX=$libInstall $ntlargs"; \
        ./configure NTL_EXCEPTIONS=on SHARED=on NTL_STD_CXX11=on NTL_SAFE_VECTORS=off TUNE=generic PREFIX=$libInstall GMP_PREFIX=$libInstall $ntlargs 2>&1 | log_install_output; \
        log_command "make"; \
        make 2>&1 | log_install_output; \
        if [ "$check" = "True" ]; then
            log_test "make check"; \
            make check 2>&1 | log_install_output; \
        fi
        log_install "make install"; \
        make install 2>&1 | log_install_output; \
        log_command "cd ../.."; \
        cd ../.. \
    )
    if [ "$remove" = "True" ]; then
        log_command "rm ntl-$ntlVers.tar.gz"
        rm ntl-$ntlVers.tar.gz || exit
    fi
    log_command "rm -rf ntl-$ntlVers"
    rm -rf ntl-$ntlVers || exit
    log_install "✅ NTL v$ntlVers instalado correctamente"
else
    log_install "⏭️  NTL ya está instalado, saltando..."
fi

# Instalación de HEXL
if [ ! -f "$libInstall/lib/libhexl.a" ]; then
    log_install "Iniciando instalación de HEXL v$hexlBranch"
    if [ ! -d "$hexlRuta" ]; then
        log_command "git clone https://github.com/intel/hexl --branch $hexlBranch $hexlRuta"
        git clone https://github.com/intel/hexl --branch $hexlBranch $hexlRuta 2>&1 | log_install_output || exit
    fi
    log_command "cd $hexlRuta"
    cd $hexlRuta || exit
    log_command "cmake -S . -B build -DCMAKE_INSTALL_PREFIX=$libInstall"
    cmake -S . -B build -DCMAKE_INSTALL_PREFIX=$libInstall 2>&1 | log_install_output || exit
    log_command "cmake --build build -j $(nproc)"
    cmake --build build -j $(nproc) 2>&1 | log_install_output || exit
    log_command "cmake --install build --prefix $libInstall"
    cmake --install build --prefix $libInstall 2>&1 | log_install_output || exit
    cd ../..
    if [ "$remove" = "True" ]; then
        log_command "rm -rf $hexlRuta"
        rm -rf $hexlRuta || exit
    else
        log_command "rm -rf $hexlRuta/hexl/build"
        rm -rf $hexlRuta/hexl/build || exit
    fi

    log_install "✅ HEXL v$hexlBranch instalado correctamente"
else
    log_install "⏭️  HEXL ya está instalado, saltando..."
fi

# Instalación de OpenFHE
log_install "Iniciando instalación de OpenFHE"
if [ ! -d "$openFHERuta" ]; then
    log_command "git clone https://github.com/openfheorg/openfhe-configurator $openFHERuta"
    git clone https://github.com/openfheorg/openfhe-configurator $openFHERuta 2>&1 | log_install_output || exit
    log_command "cd $openFHERuta"
    cd $openFHERuta || exit
else
    log_command "cd $openFHERuta && git pull"
    cd $openFHERuta && git pull 2>&1 | log_install_output || exit
fi

# Instalación de HElib
log_install "Iniciando instalación de HElib"
if [ ! -d "$helibRuta" ]; then
    log_command "git clone https://github.com/homenc/HElib $helibRuta"
    git clone https://github.com/homenc/HElib $helibRuta 2>&1 | log_install_output || exit
    log_command "cd $helibRuta"
    cd $helibRuta || exit
else
    log_command "cd $helibRuta && git pull"
    cd $helibRuta && git pull 2>&1 | log_install_output || exit
fi

if [ "$openmpInAlg" = "ON" ]; then
    algParallelism="MultiThreadFor($(nproc))"
else 
    algParallelism="SingleThread"
fi

for openFHEArg in "${openFHEArgs[@]}"; do
    log_install "Compilando OpenFHE con configuración: $openFHEArg"
    log_command "cd $openFHERuta"
    cd $openFHERuta || exit
    log_command "rm -rf ./openfhe-staging"
    rm -rf ./openfhe-staging || exit
    OPENFHE_INSTALL_DIR=$libInstall
    log_command "OPENFHE_INSTALL_DIR=$OPENFHE_INSTALL_DIR $openFHERuta/scripts/stage-openfhe-development-hexl.sh"
    OPENFHE_INSTALL_DIR=$OPENFHE_INSTALL_DIR $openFHERuta/scripts/stage-openfhe-development-hexl.sh 2>&1 | log_install_output || exit

    CMAKE_FLAGS="-DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=$libInstall -DBUILD_UNITTESTS=OFF -DBUILD_EXAMPLES=OFF -DBUILD_BENCHMARKS=OFF -DBUILD_EXTRAS=OFF -DWITH_TCM=ON -DWITH_NATIVEOPT=ON -DWITH_NTL=ON -DNTL_INCLUDE_DIR=$libInstall/include -DNTL_LIBRARIES=$libInstall/lib -DMATHBACKEND=6 $openFHEArg"
    
    log_install "Building openfhe-development with CC=$CC CXX=$CXX CMAKE_FLAGS=$CMAKE_FLAGS"

    log_command "cd $openFHERuta/openfhe-staging/openfhe-development"
    cd $openFHERuta/openfhe-staging/openfhe-development || exit
    log_command "rm -rf build"
    rm -rf build
    log_command "cmake -S . -B build $CMAKE_FLAGS"
    cmake -S . -B build $CMAKE_FLAGS 2>&1 | log_install_output || exit
    if [[ $CMAKE_FLAGS == *"-DWITH_TCM=ON"* ]]; then
        log_command "cd build"
        cd build || exit
        log_command "make tcm"
        make tcm 2>&1 | log_install_output || exit
        log_command "cd .."
        cd ..
    fi
    log_command "cmake --build build -j $(nproc)"
    cmake --build build -j $(nproc) 2>&1 | log_install_output || exit
    log_command "cmake --install build"
    cmake --install build 2>&1 | log_install_output || exit
    log_install "✅ OpenFHE compilado correctamente con $openFHEArg"

    # Crear y compilar tests de OpenFHE
    log_install "Iniciando compilación de OpenFHE_Tests para openFHEArg=$openFHEArg"
    log_command "mkdir -p $ruta_job/OpenFHE_Tests/$openFHEArg"
    mkdir -p $ruta_job/OpenFHE_Tests/$openFHEArg
    log_command "cmake -S $raizProyecto/Tests/OpenFHE_Tests -B $ruta_job/OpenFHE_Tests/$openFHEArg -DOPENFHE_ROOT=$libInstall -DENABLE_OPENMP=$openmpInAlg"
    cmake -S $raizProyecto/Tests/OpenFHE_Tests -B $ruta_job/OpenFHE_Tests/$openFHEArg -DOPENFHE_ROOT=$libInstall -DENABLE_OPENMP=$openmpInAlg 2>&1 | log_install_output || exit
    log_command "cmake --build $ruta_job/OpenFHE_Tests/$openFHEArg -j $(nproc)"
    cmake --build $ruta_job/OpenFHE_Tests/$openFHEArg -j $(nproc) 2>&1 | log_install_output || exit
    log_install "✅ OpenFHE_Tests compilado correctamente con $openFHEArg"
    exe="$ruta_job/OpenFHE_Tests/$openFHEArg/matMultTest"
    
    ejec_OpenFHE_Test() {
        local rei="$1"
        local alg="$2"
        local tamM="$3"
        local tamK="$4"
        local tamN="$5"
        log_test "[PLAN] reIt=$rei alg=$alg tamM=$tamM tamK=$tamK tamN=$tamN"
        # intentar combinaciones de ringDim/depth/otras dimensiones hasta éxito o timeout
        found_success=0
        for rdm in "${OpenFHERingDim[@]}"; do
            if [[ $found_success -ge 1 ]]; then break; fi
            log_test "[PLAN] reIt=$rei alg=$alg tamM=$tamM tamK=$tamK tamN=$tamN rdm=$rdm"
            ringDimFail=false
            for depth in "${OpenFHEMultDepth[@]}"; do
                if [[ $found_success -ge 1 || $ringDimFail == true ]]; then break; fi
                log_test "[PLAN] reIt=$rei alg=$alg tamM=$tamM tamK=$tamK tamN=$tamN rdm=$rdm depth=$depth"
                for rep in $(seq 1 $repeticiones); do
                    cmd=("$exe" --mat1 "$ficMat1" --mat2 "$ficMat2" --res "$ruta_job/OpenFHE_$ficResul" --alg "$alg" --m "$tamM" --k "$tamK" --n "$tamN" --ringDim "$rdm" --multDepth "$depth" --reIter "$rei" --rep "$rep")
                    log_command "${cmd[*]}"
                    # Ejecutar con timeout; pipe a logger y obtener el exit code del timeout
                    timeout ${timeoutSec}s "${cmd[@]}" 2>&1 | log_test_output
                    rc=${PIPESTATUS[0]}
                    case $rc in
                        0)
                            log_result "${GREEN}✅[OK] alg=$alg reIt=$rei m=$tamM k=$tamK n=$tamN rdm=$rdm depth=$depth rep=$rep${RESET}"
                            found_success=$((found_success+1))
                            ;;
                        1)
                            log_result "${YELLOW}⚠️[WARNING][BIG ERROR] alg=$alg reIt=$rei m=$tamM k=$tamK n=$tamN rdm=$rdm depth=$depth rep=$rep${RESET}"
                            ;;
                        2)
                            # Error de profundidad multiplicativa insuficiente: saltar al siguiente depth
                            log_result "${RED}❌[FAIL][DEPTH] alg=$alg reIt=$rei m=$tamM k=$tamK n=$tamN rdm=$rdm depth=$depth rep=$rep -- multiplicative depth insuficiente, probando next depth${RESET}"
                            break
                            ;;
                        3)
                            log_result "${RED}❌[FAIL][RING DIM] alg=$alg reIt=$rei m=$tamM k=$tamK n=$tamN rdm=$rdm depth=$depth rep=$rep -- probando siguiente RingDim${RESET}"
                            ringDimFail=true
                            break
                            ;;
                        124)
                            echo "${algParallelism},${tamM},${tamK},${tamN},${alg},${rei},${rep},${rdm},,${depth},,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,\"Timeout:${rc}\"" >> $ruta_job/OpenFHE_$ficResul
                            log_result "${RED}⏰[TIMEOUT] alg=$alg reIt=$rei m=$tamM k=$tamK n=$tamN rdm=$rdm depth=$depth rep=$rep (>$timeoutSec s)${RESET}"
                            return 124
                            ;;
                        134)
                            echo "${algParallelism},${tamM},${tamK},${tamN},${alg},${rei},${rep},${rdm},,${depth},,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,\"Abort:${rc}\"" >> $ruta_job/OpenFHE_$ficResul
                            log_result "${RED}❌[FAIL][ABORT] alg=$alg reIt=$rei m=$tamM k=$tamK n=$tamN rdm=$rdm depth=$depth rep=$rep -- Abort${RESET}"
                            break
                            ;;
                        137)
                            echo "${algParallelism},${tamM},${tamK},${tamN},${alg},${rei},${rep},${rdm},,${depth},,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,\"OutOfMemory:${rc}\"" >> $ruta_job/OpenFHE_$ficResul
                            log_result "${RED}❌[FAIL][OUT OF MEMORY] alg=$alg reIt=$rei m=$tamM k=$tamK n=$tamN rdm=$rdm depth=$depth rep=$rep -- Out of Memory${RESET}"
                            return 137
                            ;;
                        *)
                            echo "${algParallelism},${tamM},${tamK},${tamN},${alg},${rei},${rep},${rdm},,${depth},,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,\"Signal:${rc}\"" >> $ruta_job/HElib_$ficResul
                            log_result "${RED}❌[FAIL][$rc] alg=$alg reIt=$rei m=$tamM k=$tamK n=$tamN rdm=$rdm depth=$depth rep=$rep -- probando siguiente combinación${RESET}"
                            break
                            ;;
                    esac
                done
                if [[ $found_success -ge 1 ]]; then return 0; fi
            done
        done
        return 1
    }

    for rei in "${reiteraciones[@]}"; do
        log_test "[PLAN] reIt=$rei"
        for alg in "${OpenFHEAlg[@]}"; do
            log_test "[PLAN] reIt=$rei alg=$alg"
            for tamM in "${tamsMatriz[@]}"; do
                tamError=""
                is_rectangular=false
                for p in "${algRectangulares[@]}"; do
                    if [[ "$alg" == "$p"* ]]; then
                        is_rectangular=true
                        break
                    fi
                done

                if $is_rectangular; then
                    for tamK in "${tamsMatriz[@]}"; do
                        tamError=""
                        if [[ $rei -gt 1 ]]; then
                            ejec_OpenFHE_Test "$rei" "$alg" "$tamM" "$tamK" "$tamK"
                            if [[ $? -ne 0 ]]; then
                                if [[ $tamM -gt $tamK ]]; then
                                    tamError="m"
                                    break
                                else
                                    break
                                fi
                            fi
                        else
                            for tamN in "${tamsMatriz[@]}"; do
                                ejec_OpenFHE_Test "$rei" "$alg" "$tamM" "$tamK" "$tamN"
                                if [[ $? -ne 0 ]]; then
                                    if [[ $tamM -ge $tamK && $tamM -ge $tamN ]]; then
                                        tamError="m"
                                        break
                                    elif [[ $tamK -ge $tamM && $tamK -ge $tamN ]]; then
                                        tamError="k"
                                        break
                                    else
                                        break
                                    fi
                                fi
                            done
                        fi
                        if [[ "$tamError" == "k" ]]; then
                            break
                        fi
                    done
                else
                    ejec_OpenFHE_Test "$rei" "$alg" "$tamM" "$tamM" "$tamM"
                    if [[ $? -ne 0 ]]; then
                        break
                    fi
                fi
                if [[ "$tamError" == "m" ]]; then
                    break
                fi
            done
        done
    done
done

for helibArg in "${helibArgs[@]}"; do
    log_install "Compilando HElib con configuración: $helibArg"
    log_command "cd $helibRuta"
    cd $helibRuta || exit
    log_command "cmake -S . -B build/ -DCMAKE_INSTALL_PREFIX=$libInstall $helibArg"
    cmake -S . -B build/ -DCMAKE_INSTALL_PREFIX=$libInstall -DFETCH_GMP=OFF -DGMP_DIR=$libInstall -DNTL_DIR=$libInstall -DHEXL_DIR=$libInstall -DUSE_INTEL_HEXL=ON -DPACKAGE_BUILD=OFF -DBUILD_SHARED=OFF -DENABLE_THREADS=ON -DENABLE_TEST=ON $helibArg 2>&1 | log_install_output
    log_command "cmake --build build -j $(nproc)"
    cmake --build build -j $(nproc) 2>&1 | log_install_output
    if [ "$check" = "True" ]; then
        log_test "ctest --test-dir build -j $(nproc) --output-on-failure"
        ctest --test-dir build -j $(nproc) --output-on-failure 2>&1 | log_install_output
    fi
    log_command "cmake --install build"
    cmake --install build 2>&1 | log_install_output
    log_install "✅ HElib compilado correctamente con $helibArg"

    # Crear y compilar tests de HElib
    log_install "Iniciando compilación de HElib_Tests para HElibArg=$helibArg"
    log_command "mkdir -p $ruta_job/HElib_Tests/$helibArg"
    mkdir -p $ruta_job/HElib_Tests/$helibArg
    log_command "cmake -S $raizProyecto/Tests/HElib_Tests -B $ruta_job/HElib_Tests/$helibArg -DHELIB_ROOT=$libInstall -DENABLE_OPENMP=$openmpInAlg"
    cmake -S $raizProyecto/Tests/HElib_Tests -B $ruta_job/HElib_Tests/$helibArg -DHELIB_ROOT=$libInstall -DENABLE_OPENMP=$openmpInAlg 2>&1 | log_install_output || exit
    log_command "cmake --build $ruta_job/HElib_Tests/$helibArg -j $(nproc)"
    cmake --build $ruta_job/HElib_Tests/$helibArg -j $(nproc) 2>&1 | log_install_output || exit
    log_install "✅ HElib_Tests compilado correctamente con $helibArg"
    exe="$ruta_job/HElib_Tests/$helibArg/matMultTest"

    ejec_HElib_Test() {
        local rei="$1"
        local alg="$2"
        local tamM="$3"
        local tamK="$4"
        local tamN="$5"
        log_test "[PLAN] reIt=$rei alg=$alg tamM=$tamM tamK=$tamK tamN=$tamN"
        # intentar combinaciones de polyModulus/depth/otras dimensiones hasta éxito o timeout
        found_success=0
        for params in "${HElibParams[@]}"; do
            if [[ $found_success -ge 1 ]]; then break; fi
            # params puede ser "16384 119 2" y hay que dividir por espacios
            read -r param_m bits c <<< "$params"
            log_test "[PLAN] reIt=$rei alg=$alg tamM=$tamM tamK=$tamK tamN=$tamN param_m=$param_m bits=$bits c=$c"
            ringDimFail=false
            for precision in "${HElibPrecisions[@]}"; do
                if [[ $found_success -ge 1 || $ringDimFail == true ]]; then break; fi
                log_test "[PLAN] reIt=$rei alg=$alg tamM=$tamM tamK=$tamK tamN=$tamN param_m=$param_m bits=$bits c=$c precision=$precision"
                for rep in $(seq 1 $repeticiones); do
                    cmd=("$exe" --mat1 "$ficMat1" --mat2 "$ficMat2" --res "$ruta_job/HElib_$ficResul" --alg "$alg" --m "$tamM" --k "$tamK" --n "$tamN" --param_m "$param_m" --bits "$bits" --c "$c" --precision "$precision" --reIter "$rei" --rep "$rep")
                    log_command "${cmd[*]}"
                    # Ejecutar con timeout; pipe a logger y obtener el exit code del timeout
                    timeout ${timeoutSec}s "${cmd[@]}" 2>&1 | log_test_output
                    rc=${PIPESTATUS[0]}
                    case $rc in
                        0)
                            log_result "${GREEN}✅[OK] alg=$alg reIt=$rei m=$tamM k=$tamK n=$tamN param_m=$param_m bits=$bits c=$c precision=$precision rep=$rep${RESET}"
                            found_success=$((found_success+1))
                            ;;
                        1)
                            log_result "${YELLOW}⚠️[WARNING][BIG ERROR] alg=$alg reIt=$rei m=$tamM k=$tamK n=$tamN rdm=$rdm depth=$depth rep=$rep${RESET}"
                            break
                            ;;
                        2)
                            # Error de profundidad multiplicativa insuficiente: saltar al siguiente depth
                            log_result "${RED}❌[FAIL][DEPTH] alg=$alg reIt=$rei m=$tamM k=$tamK n=$tamN rdm=$rdm depth=$depth rep=$rep -- multiplicative depth insuficiente, probando next depth${RESET}"
                            break
                            ;;
                        3)
                            log_result "${RED}❌[FAIL][RING DIM] alg=$alg reIt=$rei m=$tamM k=$tamK n=$tamN rdm=$rdm depth=$depth rep=$rep -- probando siguiente RingDim${RESET}"
                            ringDimFail=true
                            break
                            ;;
                        124)
                            echo "${algParallelism},${tamM},${tamK},${tamN},${alg},${rei},${rep},${param_m},${bits},${c},${precision},,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,\"Timeout:${rc}\"" >> $ruta_job/HElib_$ficResul
                            log_result "${RED}⏰[TIMEOUT] alg=$alg reIt=$rei m=$tamM k=$tamK n=$tamN rdm=$rdm depth=$depth rep=$rep (>$timeoutSec s)${RESET}"
                            return 124
                            ;;
                        134)
                            echo "${algParallelism},${tamM},${tamK},${tamN},${alg},${rei},${rep},${param_m},${bits},${c},${precision},,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,\"Abort:${rc}\"" >> $ruta_job/HElib_$ficResul
                            log_result "${RED}❌[FAIL][ABORT] alg=$alg reIt=$rei m=$tamM k=$tamK n=$tamN rdm=$rdm depth=$depth rep=$rep -- Out of Memory${RESET}"
                            break
                            ;;
                        137)
                            echo "${algParallelism},${tamM},${tamK},${tamN},${alg},${rei},${rep},${param_m},${bits},${c},${precision},,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,\"OutOfMemory:${rc}\"" >> $ruta_job/HElib_$ficResul
                            log_result "${RED}❌[FAIL][OUT OF MEMORY] alg=$alg reIt=$rei m=$tamM k=$tamK n=$tamN rdm=$rdm depth=$depth rep=$rep -- Out of Memory${RESET}"
                            return 137
                            ;;
                        *)
                            echo "${algParallelism},${tamM},${tamK},${tamN},${alg},${rei},${rep},${param_m},${bits},${c},${precision},,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,\"Signal:${rc}\"" >> $ruta_job/HElib_$ficResul
                            log_result "${RED}❌[FAIL][$rc] alg=$alg reIt=$rei m=$tamM k=$tamK n=$tamN rdm=$rdm depth=$depth rep=$rep -- probando siguiente combinación${RESET}"
                            break
                            ;;
                    esac
                done
                if [[ $found_success -ge 1 ]]; then return 0; fi
            done
        done
        return 1
    }

    for rei in "${reiteraciones[@]}"; do
        log_test "[PLAN] reIt=$rei"
        for alg in "${OpenFHEAlg[@]}"; do
            log_test "[PLAN] reIt=$rei alg=$alg"
            for tamM in "${tamsMatriz[@]}"; do
                tamError=""
                is_rectangular=false
                for p in "${algRectangulares[@]}"; do
                    if [[ "$alg" == "$p"* ]]; then
                        is_rectangular=true
                        break
                    fi
                done

                if $is_rectangular; then
                    for tamK in "${tamsMatriz[@]}"; do
                        tamError=""
                        if [[ $rei -gt 1 ]]; then
                            ejec_HElib_Test "$rei" "$alg" "$tamM" "$tamK" "$tamK"
                            if [[ $? -ne 0 ]]; then
                                if [[ $tamM -gt $tamK ]]; then
                                    tamError="m"
                                    break
                                else
                                    break
                                fi
                            fi
                        else
                            for tamN in "${tamsMatriz[@]}"; do
                                ejec_HElib_Test "$rei" "$alg" "$tamM" "$tamK" "$tamN"
                                if [[ $? -ne 0 ]]; then
                                    if [[ $tamM -ge $tamK && $tamM -ge $tamN ]]; then
                                        tamError="m"
                                        break
                                    elif [[ $tamK -ge $tamM && $tamK -ge $tamN ]]; then
                                        tamError="k"
                                        break
                                    else
                                        break
                                    fi
                                fi
                            done
                        fi
                        if [[ "$tamError" == "k" ]]; then
                            break
                        fi
                    done
                else
                    ejec_HElib_Test "$rei" "$alg" "$tamM" "$tamM" "$tamM"
                    if [[ $? -ne 0 ]]; then
                        break
                    fi
                fi
                if [[ "$tamError" == "m" ]]; then
                    break
                fi
            done
        done
    done
done

exit 0