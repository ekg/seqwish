# Publishing seqwish to crates.io - Step-by-Step Guide

## Current Status

Your CLI is already complete and uses `clap` for argument parsing. The command is:

```bash
seqwish -s sequences.fa -p alignments.paf -g output.gfa [OPTIONS]
```

## Main Blocking Issue: iitree-rs Dependency

**Current problem**: seqwish depends on iitree-rs via git, but crates.io requires all dependencies to be published.

```toml
# Current (blocks publication):
iitree-rs = { git = "https://github.com/pangenome/iitree-rs", rev = "41cb916..." }

# Required for publication:
iitree-rs = "0.1.0"
```

## Two-Step Publication Process

### Step 1: Publish iitree-rs First

iitree-rs is READY TO PUBLISH! It has all required metadata in Cargo.toml.

**From `/home/erik/iitree-rs` directory:**

```bash
cd /home/erik/iitree-rs

# 1. Get your crates.io API token
#    Go to https://crates.io/me (must be logged in)
#    Click "Account Settings" → "API Tokens" → "New Token"
#    Give it a name like "iitree-rs-publish"
#    Copy the token

# 2. Save your token (one-time setup)
cargo login
# Paste your token when prompted
# Token will be saved to ~/.cargo/credentials.toml

# 3. Test the publication (doesn't actually upload)
cargo publish --dry-run

# 4. If dry-run succeeds, publish for real
cargo publish

# Note: You can only publish each version ONCE - no going back!
```

### Step 2: Update seqwish and Publish

**From `/home/erik/seqwish` directory:**

```bash
cd /home/erik/seqwish

# 1. Update Cargo.toml to use published iitree-rs
#    Change this line:
#    iitree-rs = { git = "https://github.com/pangenome/iitree-rs", rev = "41cb916..." }
#    To this:
#    iitree-rs = "0.1.0"

# 2. Update lockfile
cargo update -p iitree-rs

# 3. Test that everything still builds
cargo build --release
cargo test

# 4. Test publication (doesn't upload)
cargo publish --dry-run

# 5. If successful, publish
cargo publish
```

## Pre-Flight Checklist

### For iitree-rs:
- ✅ Has Cargo.toml metadata (name, version, authors, description, license, repository)
- ✅ Has keywords and categories
- ✅ All dependencies are from crates.io
- ✅ README.md exists
- ⚠️  Need to check if LICENSE file exists
- ⚠️  Need to verify it builds cleanly

### For seqwish:
- ✅ Has Cargo.toml metadata
- ✅ Has comprehensive README.md
- ✅ Has API documentation in lib.rs
- ❌ Needs LICENSE file
- ❌ Depends on unpublished iitree-rs (must fix)
- ✅ All other dependencies from crates.io

## Quick Setup Commands

```bash
# Check if you're already logged in
cat ~/.cargo/credentials.toml 2>/dev/null

# If not, get your token from https://crates.io/me
# Then:
cargo login

# Test iitree-rs
cd /home/erik/iitree-rs
cargo publish --dry-run

# If it complains about LICENSE, create one:
cat > LICENSE << 'EOF'
MIT License

Copyright (c) 2024 Erik Garrison

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
EOF
```

## What cargo publish Does

1. **Packages** your crate into a `.crate` file (compressed tarball)
2. **Verifies** by building in a clean environment
3. **Uploads** to crates.io (if not --dry-run)
4. **Indexes** so others can `cargo install seqwish`

## After Publishing

Once published:
- Anyone can install: `cargo install seqwish`
- Can be used as dependency: `seqwish = "0.1.0"`
- Version is PERMANENT - can't delete or modify
- Can publish new versions: bump version number in Cargo.toml

## Common Issues

**"cannot publish with dirty repository"**
```bash
# Commit your changes first, or use:
cargo publish --allow-dirty  # Not recommended
```

**"failed to verify package"**
```bash
# It's building in a clean environment and failing
# Usually means a file is missing from version control
cargo publish --no-verify  # Skip verification (not recommended)
```

**"crate name already taken"**
```bash
# Someone else owns that name on crates.io
# Need to rename in Cargo.toml or request transfer
```

## Versioning Strategy

After 0.1.0, follow [Semantic Versioning](https://semver.org/):

- **0.1.0** → **0.1.1**: Bug fixes (patch)
- **0.1.0** → **0.2.0**: New features, backward compatible (minor)
- **0.1.0** → **1.0.0**: Breaking changes or "production ready" (major)

Update version in Cargo.toml, then `cargo publish` again.

## Next Steps

1. Check LICENSE files for both crates
2. Do dry-run publishes to catch issues
3. Publish iitree-rs
4. Update seqwish Cargo.toml
5. Publish seqwish
6. Celebrate! 🎉

## Resources

- Crates.io Publishing Guide: https://doc.rust-lang.org/cargo/reference/publishing.html
- Your crates.io account: https://crates.io/me
- Manage published crates: https://crates.io/crates/[crate-name]/settings
