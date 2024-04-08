cpmaddpackage(
        NAME libassert
        GIT_REPOSITORY "git@github.com:jeremy-rifkin/libassert.git"
        GIT_TAG main)

sl_target_link_system_libraries(${PROJECT_NAME}
        INTERFACE
        # utils
        libassert::assert
)

# Packages for tests and examples

if (NOT PROJECT_IS_TOP_LEVEL)
    return()
endif ()

cpmaddpackage(
        NAME range-v3
        GIT_REPOSITORY "git@github.com:ericniebler/range-v3.git"
        GIT_TAG 0.12.0)
cpmaddpackage(
        NAME fmt
        GIT_REPOSITORY "git@github.com:fmtlib/fmt.git"
        GIT_TAG 10.1.1)

sl_target_link_system_libraries(${PROJECT_NAME}
        INTERFACE
        range-v3::range-v3
        fmt::fmt
)
