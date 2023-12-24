cpmaddpackage(
        NAME range-v3
        GIT_REPOSITORY "git@github.com:ericniebler/range-v3.git"
        GIT_TAG 0.12.0)


sl_add_example(${PROJECT_NAME} fft_exploration)
sl_target_link_system_libraries(${PROJECT_NAME}_fft_exploration
        PRIVATE
        range-v3::range-v3
)

