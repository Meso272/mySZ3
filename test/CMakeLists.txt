

add_executable(sz_huffman_test sz_huffman_test.cpp)
target_include_directories(sz_huffman_test PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_link_libraries(sz_huffman_test SZ3)

add_executable(sz_test sz_compress_test.cpp)
target_include_directories(sz_test PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_link_libraries(sz_test SZ3)

add_executable(sz_demo sz_demo.cpp)
target_include_directories(sz_demo PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_link_libraries(sz_demo SZ3)

add_executable(sz_backend sz_backend.cpp)
target_include_directories(sz_backend PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_link_libraries(sz_backend SZ3)

add_executable(sz_backend_ae sz_backend_ae.cpp)
target_include_directories(sz_backend_ae PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_link_libraries(sz_backend_ae SZ3)