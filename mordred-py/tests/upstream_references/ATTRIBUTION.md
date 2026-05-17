# Vendored upstream reference corpus

The files in this directory are vendored verbatim from
[mordredcommunity/mordred](https://github.com/mordredcommunity/mordred)
(`mordred/tests/references/`) and serve as the parity gate for `mordred-rs`:
expected descriptor values that the Rust port must reproduce to be considered
drop-in compatible.

## Provenance

- **Source repo:** https://github.com/mordredcommunity/mordred
- **Subdirectory:** `mordred/tests/references/`
- **Snapshot:** vendored 2026-05-13

## License

Upstream is licensed under BSD-3-Clause:

```
Copyright (c) 2023, Jackson Burns and the mordredcommunity Team
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.
3. Neither the name of the copyright holder nor the names of its contributors
   may be used to endorse or promote products derived from this software
   without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
POSSIBILITY OF SUCH DAMAGE.
```

## Refreshing

To re-vendor from upstream:

```sh
git clone https://github.com/mordredcommunity/mordred /tmp/mordred-upstream
rm -rf mordred-py/tests/upstream_references
cp -R /tmp/mordred-upstream/mordred/tests/references \
      mordred-py/tests/upstream_references
# preserve this ATTRIBUTION.md (re-create from this template)
```

## Why not symlink or git-submodule?

A snapshot in-tree:

- Keeps the parity contract stable across upstream churn (rare, but possible).
- Avoids the symlink portability issues across macOS/Linux/Windows CI runners.
- Lets contributors reproduce the test gate without a second clone.

The trade-off is manual refresh, which is acceptable because upstream mordred
is largely in maintenance mode.
