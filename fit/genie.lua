solution "fitLTC"
   configurations { "Debug", "Release" }

   project "fitLTC"
      kind "ConsoleApp"
      language "C++"
      files { "**.h", "**.cpp", "**.c" }

      includedirs {
         "../external/CImg",
         "../external/glm"
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
