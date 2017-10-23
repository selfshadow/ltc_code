solution "prefilter"
   configurations { "Debug", "Release" }

   project "prefilter"
      kind "ConsoleApp"
      language "C++"
      files { "**.h", "**.cpp", "**.c" }

      includedirs {
         "../external/CImg",
         "../external/glm",
         "../external/stb"
      }

      defines { 
         "cimg_display=0"
      }

      configuration "Debug"
         targetdir "bin"
         defines { "DEBUG" }
         flags { "Symbols" }

      configuration "Release"
         targetdir "bin"
         defines { "NDEBUG" }
         flags { "Optimize" }
