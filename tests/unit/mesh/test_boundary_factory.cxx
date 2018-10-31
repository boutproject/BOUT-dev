#include "gtest/gtest.h"

#include "boundary_factory.hxx"
#include "boundary_region.hxx"

#include "test_extras.hxx"

/// Global mesh
extern Mesh *mesh;

class TestBoundary : public BoundaryOp {
public:
  BoundaryOp *clone(BoundaryRegion *UNUSED(region), const std::list<std::string> &args,
                    const std::map<std::string, std::string> &keywords) override {
    auto *testboundary = new TestBoundary();
    testboundary->args = args;
    testboundary->keywords = keywords;
    return testboundary;
  }
  std::list<std::string> args;
  std::map<std::string, std::string> keywords;

  void apply(Field2D &UNUSED(f)) override {}
  void apply(Field3D &UNUSED(f)) override {}
};

TEST(BoundaryFactoryTests, IsSingleton) {
  BoundaryFactory* fac1 = BoundaryFactory::getInstance();
  BoundaryFactory* fac2 = BoundaryFactory::getInstance();

  EXPECT_EQ(fac1, fac2);
}

TEST(BoundaryFactoryTests, CreateTestBoundary) {
  mesh = new FakeMesh(3, 3, 3);
  
  BoundaryFactory* fac = BoundaryFactory::getInstance();
  fac->add(new TestBoundary(), "testboundary");
  
  BoundaryRegionXIn region("test_region", 0, 1, mesh);

  // Check no brackets
  
  auto *boundary = fac->create("testboundary", &region);

  EXPECT_TRUE( boundary != nullptr );
  
  EXPECT_EQ(typeid(*boundary), typeid(TestBoundary));

  delete boundary;
  
  // Positional arguments

  boundary = fac->create("testboundary(a, 1)", &region);
  EXPECT_TRUE( boundary != nullptr );

  TestBoundary *tb = dynamic_cast<TestBoundary*>(boundary);

  EXPECT_EQ( tb->args.size(), 2 );
  EXPECT_EQ( tb->args.front(), "a" );
  EXPECT_EQ( tb->args.back(), "1" );
  EXPECT_TRUE( tb->keywords.empty() );
  
  delete boundary;

  // Test keywords
  boundary = fac->create("testboundary(key=1, b=value)", &region);
  EXPECT_TRUE( boundary != nullptr );

  tb = dynamic_cast<TestBoundary*>(boundary);
  
  EXPECT_TRUE( tb->args.empty() );
  EXPECT_EQ( tb->keywords.size(), 2 );
  EXPECT_EQ( tb->keywords.at("key"), "1" );
  EXPECT_EQ( tb->keywords.at("b"), "value" );

  delete boundary;

  // Mix of positional args and keywords
  boundary = fac->create("testboundary(0.23, key =1+2 , something(),b=value ,a + sin(1.2))", &region);
  EXPECT_TRUE( boundary != nullptr );

  tb = dynamic_cast<TestBoundary*>(boundary);

  std::list<std::string> args_expected = {"0.23", "something()", "a + sin(1.2)"};
  
  EXPECT_EQ( tb->args, args_expected );
  
  EXPECT_EQ( tb->keywords.size(), 2 );
  EXPECT_EQ( tb->keywords.at("key"), "1+2" );
  EXPECT_EQ( tb->keywords.at("b"), "value" );
}
