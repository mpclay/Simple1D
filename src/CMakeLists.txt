# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
#
# Copyright (C) 2014 by Matthew Clay <mpclay@gmail.com>

# Core source files that are program independent.
SET(CORE_SOURCE_FILES
    ${CMAKE_CURRENT_SOURCE_DIR}/Parameters.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/BoundaryConditions.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/Flux.F90
    CACHE INTERNAL "Core source files.")

# Source files for the Simple1D program.
SET(SIMPLE1D_SOURCE_FILES
    ${CMAKE_CURRENT_SOURCE_DIR}/Simple1D.F90
    CACHE INTERNAL "Program file for Simple1D.")

