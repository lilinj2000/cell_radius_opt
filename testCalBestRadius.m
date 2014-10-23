function test_suite = testCalBestRadius
initTestSuite;

function testCalBestRadius_k_1

geo_radius = [ 1:5 ];
r_radius =  geo_radius;

expect_radius = r_radius;
expect_similar_best = 1;
expect_k = 1;

[radius, similar_best, k] = calBestRadius(geo_radius, r_radius);

assertVectorsAlmostEqual(radius, expect_radius);
assertEqual(similar_best, expect_similar_best);
assertElementsAlmostEqual(k, expect_k);

function testCalBestRadius_k_1_2

geo_radius = [1:10];
r_radius = geo_radius*1.2;

expect_radius = r_radius;
expect_similar_best = 1;
expect_k = 1.2;

[radius, similar_best, k] = calBestRadius(geo_radius, r_radius);

assertVectorsAlmostEqual(radius, expect_radius);
assertEqual(similar_best, expect_similar_best);
assertElementsAlmostEqual(k, expect_k);