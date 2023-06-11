
add_subdirectory(lib)

add_executable(test1 test1.cpp)
add_executable(test2 test2.cpp)
add_executable(test3 test3.cpp)
add_executable(test4 test4.cpp)

target_link_libraries(test1 ${PROJECT_NAME})
target_link_libraries(test2 ${PROJECT_NAME})
target_link_libraries(test3 ${PROJECT_NAME})
target_link_libraries(test4 ${PROJECT_NAME})

install(TARGETS
        test1
        test2
        test3
        test4
        RUNTIME DESTINATION .
        LIBRARY DESTINATION .
        ARCHIVE DESTINATION .)
