#!/usr/bin/env bash
set -e

# ------------------------------------------------------------------------------
# 1) Build your program (adjust 'release' to whatever target you want).
#    We'll assume your existing Makefile has a 'release' target that
#    produces 'build/util_exec'.
# ------------------------------------------------------------------------------
echo "==> Building the project in release mode..."
make release

# ------------------------------------------------------------------------------
# 2) Create a distribution folder: build-dist/
# ------------------------------------------------------------------------------
DIST_DIR="build-dist"
rm -rf "$DIST_DIR"
mkdir -p "$DIST_DIR/libs"

# Copy the main executable from 'build/' into 'build-dist/'
cp build/util_exec "$DIST_DIR/"

# ------------------------------------------------------------------------------
# 3) Function to copy any /usr/local or /opt/homebrew libraries used by a binary
#    into build-dist/libs, and fix references in that binary to use @executable_path/libs/...
# ------------------------------------------------------------------------------
copy_and_fix_executable() {
  local binary="$1"

  # Use otool -L to find all references to /usr/local or /opt/homebrew
  local deps
  deps=$(otool -L "$binary" | awk '/\/usr\/local|\/opt\/homebrew/ {print $1}')
  
  for dep in $deps; do
    local depname
    depname=$(basename "$dep")
    # If this library isn't already copied, do so
    if [ ! -f "$DIST_DIR/libs/$depname" ]; then
      echo "==> Copying $dep -> $DIST_DIR/libs/$depname"
      cp "$dep" "$DIST_DIR/libs/$depname"
    fi

    echo "==> install_name_tool -change $dep @executable_path/libs/$depname $binary"
    install_name_tool -change "$dep" "@executable_path/libs/$depname" "$binary"
  done
}

# ------------------------------------------------------------------------------
# 4) Function to fix references inside all .dylibs in build-dist/libs/
#    so they don't point back to /usr/local or /opt/homebrew.
#    Instead, we use @rpath/<libname>.
# ------------------------------------------------------------------------------
fix_libraries_in_libs_folder() {
  for dylib in "$DIST_DIR"/libs/*.dylib; do
    [ -f "$dylib" ] || continue  # skip if not a file

    echo "==> Processing library: $dylib"

    # 4a) Fix each /usr/local or /opt/homebrew reference
    local inner_deps
    inner_deps=$(otool -L "$dylib" | awk '/\/usr\/local|\/opt\/homebrew/ {print $1}')
    for dep in $inner_deps; do
      local depname
      depname=$(basename "$dep")

      # Copy it if we haven't yet
      if [ ! -f "$DIST_DIR/libs/$depname" ]; then
        echo "     copying $dep -> $DIST_DIR/libs/$depname"
        cp "$dep" "$DIST_DIR/libs/$depname"
      fi

      echo "     install_name_tool -change $dep @rpath/$depname $dylib"
      install_name_tool -change "$dep" "@rpath/$depname" "$dylib"
    done

    # 4b) Set this library's own ID to @rpath/<filename>, so other libs referencing it find it via rpath
    local libfile
    libfile=$(basename "$dylib")
    echo "     install_name_tool -id @rpath/$libfile $dylib"
    install_name_tool -id "@rpath/$libfile" "$dylib"
  done
}

# ------------------------------------------------------------------------------
# 5) First pass: fix references in the main executable
# ------------------------------------------------------------------------------
echo "==> First pass: fix references in main executable"
copy_and_fix_executable "$DIST_DIR/util_exec"

# ------------------------------------------------------------------------------
# 6) Fix references inside each .dylib that got copied so far
# ------------------------------------------------------------------------------
echo
echo "==> Fix references inside *.dylib in $DIST_DIR/libs"
fix_libraries_in_libs_folder

# ------------------------------------------------------------------------------
# 7) Second pass: It's often necessary to fix the main executable again if
#    we discovered new libs in step 6 (for example libgfortran.5.dylib, etc.)
# ------------------------------------------------------------------------------
echo
echo "==> Second pass: fix references in main executable again (in case new libs appeared)"
copy_and_fix_executable "$DIST_DIR/util_exec"

# ------------------------------------------------------------------------------
# 8) (Optional) Re-check references inside libs one more time to catch newly copied libs
#    We'll skip it here for brevity. If you see leftover references, you can repeat step 6.
# ------------------------------------------------------------------------------
echo
echo "==> Final pass: fix references inside libs again (optional, but safer)."
fix_libraries_in_libs_folder

# ------------------------------------------------------------------------------
# 9) Show final references
# ------------------------------------------------------------------------------
echo
echo "==> Final check: references in $DIST_DIR/util_exec"
otool -L "$DIST_DIR/util_exec"
echo
echo "==> Final check: references in $DIST_DIR/libs/*.dylib"
otool -L "$DIST_DIR"/libs/*.dylib || true

# ------------------------------------------------------------------------------
# 10) Zip it up
# ------------------------------------------------------------------------------
echo
echo "==> Creating util_exec-macos.zip..."
zip -r util_exec-macos.zip "$DIST_DIR"

echo
echo "All done! Distribute 'util_exec-macos.zip' to your users."
echo "They can unzip, 'cd build-dist', and run './util_exec'."
