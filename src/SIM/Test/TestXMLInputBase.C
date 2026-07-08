//==============================================================================
//!
//! \file TestXMLInputBase.C
//!
//! \date Nov 28 2023
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Tests for base class for XML input parsing.
//!
//==============================================================================

#include "XMLInputBase.h"

#include "Catch2Support.h"

#include <string>
#include <vector>
#include <tinyxml2.h>


namespace {

class TestXMLInput : public XMLInputBase
{
public:
  bool parse(const tinyxml2::XMLElement* elem) override
  {
    auto parseElem = [&content=strings](const tinyxml2::XMLElement* elm)
    {
      if (!elm) return false;

      content.push_back(elm->Value());
      for (const tinyxml2::XMLAttribute* att = elm->FirstAttribute();
           att; att = att->Next())
      {
        content.push_back(att->Name());
        content.push_back(att->Value());
      }
      if (elm->GetText())
        content.push_back(elm->GetText());

      return true;
    };

    if (parseElem(elem))
      for (const tinyxml2::XMLElement* child = elem->FirstChildElement();
           child; child = child->NextSiblingElement())
        parseElem(child);

    return true;
  }

  std::vector<std::string> strings;
};

}


TEST_CASE("TestXMLInputBase.IncludeFiles")
{
  const std::vector<std::string> ref = {
    "boundaryconditions",
    "dirichlet", "set", "foo", "comp", "12",
    "dirichlet", "set", "bar", "comp", "12", "type", "expression", "a*b*c",
    "neumann", "set", "foobar", "type", "constant", "1.0",
    "neumann", "set", "barbar", "type", "constant", "2.0",
    "someothertag", "is_here"
  };

  TestXMLInput x;
  REQUIRE(x.readXML("src/SIM/Test/with_include.xml"));
  REQUIRE(x.strings == ref);
}
