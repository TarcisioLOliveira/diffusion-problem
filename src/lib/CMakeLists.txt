set(SOURCES
    eigen.cpp
    mesh.cpp
    Q4.cpp
    sparse_matrix.cpp
    window.cpp
)

add_library(${PROJECT_NAME} ${SOURCES})

target_link_libraries(${PROJECT_NAME} PUBLIC ${LAPACKE_LIBRARIES} cblas ${BLAS_LIBRARIES} ${LAPACKE_LIBRARIES} Eigen3::Eigen sfml-graphics)

target_link_libraries(${PROJECT_NAME} PUBLIC OpenMP::OpenMP_CXX OpenMP::OpenMP_C OpenMP::OpenMP_Fortran)

install(TARGETS ${PROJECT_NAME} RUNTIME DESTINATION .
                                LIBRARY DESTINATION .
                                ARCHIVE DESTINATION .)
