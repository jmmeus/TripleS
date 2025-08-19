"""Test requirements.txt file."""
from pathlib import Path


def test_requirements_file_exists() -> None:
    """Test that requirements.txt exists."""
    requirements_path = Path(__file__).parent.parent / "requirements.txt"
    assert requirements_path.exists()


def test_requirements_has_expected_dependencies() -> None:
    """Test that requirements.txt contains expected dependencies."""
    requirements_path = Path(__file__).parent.parent / "requirements.txt"
    content = requirements_path.read_text().strip()
    
    dependencies = [line.strip() for line in content.split('\n') if line.strip()]
    
    expected_deps = {"matplotlib", "numpy", "pytest"}
    actual_deps = set(dependencies)
    
    assert expected_deps.issubset(actual_deps), f"Missing dependencies: {expected_deps - actual_deps}"