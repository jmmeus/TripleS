"""Test README.md contributors section."""
from pathlib import Path


def test_readme_exists() -> None:
    """Test that README.md exists."""
    readme_path = Path(__file__).parent.parent / "README.md"
    assert readme_path.exists()


def test_contributors_section_exists() -> None:
    """Test that Contributors section exists in README."""
    readme_path = Path(__file__).parent.parent / "README.md"
    content = readme_path.read_text()
    assert "## Contributors" in content


def test_contributors_list_not_empty() -> None:
    """Test that contributors list contains names."""
    readme_path = Path(__file__).parent.parent / "README.md"
    content = readme_path.read_text()
    
    lines = content.split('\n')
    contributors_section_found = False
    
    for i, line in enumerate(lines):
        if line.strip() == "## Contributors":
            contributors_section_found = True
            # Check next non-empty line contains contributors
            if i + 2 < len(lines):
                contributors_line = lines[i + 2].strip()
                assert len(contributors_line) > 0, "Contributors list is empty"
                assert "," in contributors_line, "Expected comma-separated contributors"
            break
    
    assert contributors_section_found, "Contributors section not found"