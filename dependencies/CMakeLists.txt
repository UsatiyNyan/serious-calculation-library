cpmaddpackage(
        NAME libassert
        GIT_REPOSITORY "git@github.com:jeremy-rifkin/libassert.git"
        GIT_TAG main)

sl_target_link_system_libraries(${PROJECT_NAME}
        INTERFACE
        # utils
        assert
)
