function run_tests
%RUN_TESTS Execute the package regression tests.

tests = {@test_scalar_exact,@test_stationary_limit,@test_symmetry_psd, ...
    @test_tolerance_response};
fprintf('Running %d tests...\n',numel(tests));
for k = 1:numel(tests)
    name = func2str(tests{k});
    fprintf('  %-30s',name);
    tests{k}();
    fprintf('PASS\n');
end
fprintf('All tests passed.\n');
end
