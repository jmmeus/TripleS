import pygame
from pygame.locals import *
from OpenGL.GL import *
from OpenGL.GLU import *
import math
import numpy as np
from typing import List, Tuple, Optional, Union, Any, Sequence

class Camera3D:
    '''
    Camera 3D model.
    Handles camera rotation controls and camera velocity.
    '''

    def __init__(self, position: Optional[List[float]] = None, target: Optional[List[float]] = None, up: Optional[List[float]] = None) -> None:
        # Input validation
        if position is not None:
            if not isinstance(position, list) or len(position) != 3 or not all(isinstance(x, (int, float)) for x in position):
                raise ValueError("Position must be a list of 3 numeric values")
        if target is not None:
            if not isinstance(target, list) or len(target) != 3 or not all(isinstance(x, (int, float)) for x in target):
                raise ValueError("Target must be a list of 3 numeric values")
        if up is not None:
            if not isinstance(up, list) or len(up) != 3 or not all(isinstance(x, (int, float)) for x in up):
                raise ValueError("Up must be a list of 3 numeric values")
            
        if position is None:
            position = [0, 0, 0]
        if target is None:
            target = [0, 0, -1]
        if up is None:
            up = [0, 1, 0]
            
        self.position: np.ndarray = np.array(position, dtype=float)
        self.target: np.ndarray = np.array(target, dtype=float)
        self.up: np.ndarray = np.array(up, dtype=float)

        # Camera rotation angles - start looking straight down -Z axis
        self.yaw: float = 0.0  # 0 degrees to look down -Z axis
        self.pitch: float = 0.0  # Level

        # Movement speeds
        self.movement_speed: float = 5.0
        self.mouse_sensitivity: float = 0.002
        self.scroll_speed: float = 2.0

        # Mouse state
        self.last_mouse_pos: Tuple[int, int] = (0, 0)
        self.mouse_captured: bool = False

        # Calculate initial forward vector
        self.update_vectors()

    def update_vectors(self) -> None:
        """Update camera vectors based on yaw and pitch"""
        # Calculate forward vector from yaw and pitch
        # Adjusted so yaw=0 looks down -Z axis
        self.forward: np.ndarray = np.array([
            -math.sin(self.yaw) * math.cos(self.pitch),
            math.sin(self.pitch),
            -math.cos(self.yaw) * math.cos(self.pitch)
        ])

        # Calculate right vector (cross product of forward and world up)
        world_up: np.ndarray = np.array([0, 1, 0])
        self.right: np.ndarray = np.cross(self.forward, world_up)
        if np.linalg.norm(self.right) > 0:
            self.right = self.right / np.linalg.norm(self.right)

        # Calculate up vector
        self.up = np.cross(self.right, self.forward)
        if np.linalg.norm(self.up) > 0:
            self.up = self.up / np.linalg.norm(self.up)

    def handle_mouse_motion(self, mouse_pos: Tuple[int, int], mouse_buttons: Tuple[bool, bool, bool]) -> None:
        """Handle mouse movement for camera rotation"""
        if not isinstance(mouse_pos, tuple) or len(mouse_pos) != 2 or not all(isinstance(x, int) for x in mouse_pos):
            raise ValueError("mouse_pos must be a tuple of 2 integers")
        if not isinstance(mouse_buttons, tuple) or len(mouse_buttons) != 3 or not all(isinstance(x, bool) for x in mouse_buttons):
            raise ValueError("mouse_buttons must be a tuple of 3 boolean values")
        if not self.mouse_captured:
            self.last_mouse_pos = mouse_pos
            return

        # Get relative mouse movement
        mouse_rel = pygame.mouse.get_rel()
        dx = mouse_rel[0]
        dy = mouse_rel[1]

        # Update rotation angles
        self.yaw -= dx * self.mouse_sensitivity
        self.pitch -= dy * self.mouse_sensitivity  # Negative for natural feel

        # Clamp pitch to prevent camera flip
        self.pitch = max(-math.pi / 2 + 0.1, min(math.pi / 2 - 0.1, self.pitch))

        # Update camera vectors
        self.update_vectors()

        self.last_mouse_pos = mouse_pos

    def handle_mouse_wheel(self, wheel_y: float) -> None:
        """Handle mouse wheel for zoom/forward movement"""
        if not isinstance(wheel_y, (int, float)):
            raise ValueError("wheel_y must be a numeric value")
        self.position += self.forward * wheel_y * self.scroll_speed

    def handle_keyboard(self, keys: Sequence[bool], dt: float) -> None:
        """Handle keyboard input for movement"""
        if not hasattr(keys, '__getitem__'):
            raise ValueError("keys must be a sequence that supports indexing")
        if not isinstance(dt, (int, float)) or dt < 0:
            raise ValueError("dt must be a non-negative numeric value")
        velocity = self.movement_speed * dt

        # WASD movement
        if keys[pygame.K_w]:
            self.position += self.forward * velocity
        if keys[pygame.K_s]:
            self.position -= self.forward * velocity
        if keys[pygame.K_a]:
            self.position -= self.right * velocity
        if keys[pygame.K_d]:
            self.position += self.right * velocity

        # QE for up/down movement
        if keys[pygame.K_q]:
            self.position += self.up * velocity
        if keys[pygame.K_e]:
            self.position -= self.up * velocity

    def toggle_mouse_capture(self) -> None:
        """Toggle mouse capture for camera control"""
        self.mouse_captured = not self.mouse_captured
        if self.mouse_captured:
            pygame.mouse.set_visible(False)
            pygame.event.set_grab(True)
            pygame.mouse.get_rel()  # Clear any pending mouse movement
        else:
            pygame.mouse.set_visible(True)
            pygame.event.set_grab(False)

    def apply_camera_transform(self) -> None:
        """Apply camera transformation to OpenGL"""
        # Calculate look-at point
        look_at = self.position + self.forward

        # Use gluLookAt to set up the camera
        gluLookAt(
            self.position[0], self.position[1], self.position[2],  # Eye position
            look_at[0], look_at[1], look_at[2],  # Look-at point
            self.up[0], self.up[1], self.up[2]  # Up vector
        )


class Renderer:
    def __init__(self) -> None:
        pygame.init()
        self.fullscreen = False
        # If fullscreen then use correct fullscreen openGL settings
        if self.fullscreen:
            self.display = pygame.display.list_modes()[0]
            self.screen = pygame.display.set_mode(self.display, DOUBLEBUF | OPENGL | FULLSCREEN)
        # Else use 1200p resolution. (Can include more resolution settings in future).
        else:
            self.display = (1200, 800)
            self.screen = pygame.display.set_mode(self.display, DOUBLEBUF | OPENGL)

        pygame.display.set_caption("3D Camera Controls")
        self.clock = pygame.time.Clock()
        self.width = self.display[0]
        self.height = self.display[1]

        # Initialize font for text overlay
        pygame.font.init()
        self.font = pygame.font.Font(None, 24)

        # Initialize OpenGL
        self.init_opengl()

        # Initialize camera - start at (0, 0, 5) looking down -Z axis
        self.camera = Camera3D(position=[0, 0, 5])

        # Sample 3D points to render (a simple cube)
        self.cube_vertices = [
            [-1, -1, -1], [1, -1, -1], [1, 1, -1], [-1, 1, -1],  # Back face
            [-1, -1, 1], [1, -1, 1], [1, 1, 1], [-1, 1, 1]  # Front face
        ]

        # Cube edges for wireframe rendering
        self.cube_edges = [
            (0, 1), (1, 2), (2, 3), (3, 0),  # Back face
            (4, 5), (5, 6), (6, 7), (7, 4),  # Front face
            (0, 4), (1, 5), (2, 6), (3, 7)  # Connecting edges
        ]

    def init_opengl(self) -> None:
        """Initialize OpenGL settings"""
        # Set clear color to dark blue
        glClearColor(0.0, 0.1, 0.2, 1.0)

        # Enable depth testing (ensures objects closer to the camera are drawn in front of further objects)
        glEnable(GL_DEPTH_TEST)

        # Set up lighting
        glEnable(GL_LIGHTING)
        glEnable(GL_LIGHT0) # Enable first light source (OpenGl supports 8 light sources up to LIGHT7)
        glEnable(GL_COLOR_MATERIAL) # Allow object colors to interact with lighting
        # GL_FRONT_AND_BACK applies to both front and back faces of polygons GL_AMBIENT_AND_DIFFUSE means object
        # colors affect both ambient and diffuse lighting
        glColorMaterial(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE)

        # Light properties
        light_position = [5.0, 5.0, 5.0, 1.0] # Where the light is located
        light_ambient = [0.2, 0.2, 0.2, 1.0] # Background lighting (like moonlight)
        light_diffuse = [0.8, 0.8, 0.8, 1.0] # Main directional lighting (like sunlight)
        light_specular = [1.0, 1.0, 1.0, 1.0] # Shiny reflections (like highlights on metal)

        glLightfv(GL_LIGHT0, GL_POSITION, light_position)
        glLightfv(GL_LIGHT0, GL_AMBIENT, light_ambient)
        glLightfv(GL_LIGHT0, GL_DIFFUSE, light_diffuse)
        glLightfv(GL_LIGHT0, GL_SPECULAR, light_specular)

        # Set up perspective projection
        glMatrixMode(GL_PROJECTION) #GL_PROJECTION: Controls how 3D coordinates are projected to 2D screen
        glLoadIdentity()
        # params for gluPerspective: Field of view in degrees (wider = more fish-eye effect),
        # aspect_ratio: Width/height ratio to prevent stretching
        # Near clipping plane: (objects closer than this aren't drawn)
        # Far clipping plane: (objects farther than this aren't drawn)
        gluPerspective(45, (self.display[0] / self.display[1]), 0.1, 100.0)

        # Switch back to modelview matrix
        glMatrixMode(GL_MODELVIEW) # GL_MODELVIEW: Controls camera position and object transformations
    def draw_planet_center(self, x: float, y: float, z: float) -> None:
        """Draw a point at planet center"""
        if not all(isinstance(coord, (int, float)) for coord in [x, y, z]):
            raise ValueError("Coordinates x, y, z must be numeric values")
        glDisable(GL_LIGHTING)
        glColor3f(1.0, 1.0, 0.0)  # Yellow color
        glPointSize(10)
        glBegin(GL_POINTS)
        glVertex3f(x, y, z)
        glEnd()
        glEnable(GL_LIGHTING)

    def draw_orbit(self, orbit_radius: float) -> None:
        """Draw orbit circle"""
        if not isinstance(orbit_radius, (int, float)) or orbit_radius <= 0:
            raise ValueError("orbit_radius must be a positive numeric value")
        glDisable(GL_LIGHTING)
        glColor3f(0.5, 0.5, 1.0)  # Light blue color for orbit
        glBegin(GL_LINE_LOOP)
        for i in range(360):
            angle = math.radians(i)
            x = orbit_radius * math.cos(angle)
            z = orbit_radius * math.sin(angle)
            glVertex3f(x, 0, z)
        glEnd()
        glEnable(GL_LIGHTING)

    def generate_solar_system(self, planets: Optional[List[Any]] = None) -> None:
        """Generate solar system objects"""
        if planets is None:
            planets = []
        if not isinstance(planets, list):
            raise ValueError("planets must be a list")
        # Draw sun at center
        glPushMatrix()
        glColor3f(1.0, 1.0, 0.0)  # Yellow sun
        gluSphere(gluNewQuadric(), 0.8, 20, 20)
        self.draw_planet_center(0, 0, 0)
        glPopMatrix()

        # Draw planets at different orbital distances
        planet_data = [
            (3.0, [0.8, 0.4, 0.2], 0.3),  # Mercury-like
            (5.0, [1.0, 0.7, 0.3], 0.4),  # Venus-like
            (7.0, [0.3, 0.5, 1.0], 0.4),  # Earth-like
            (9.0, [1.0, 0.3, 0.3], 0.35),  # Mars-like
        ]

        for orbit_radius, color, planet_size in planet_data:
            glPushMatrix()

            # Draw orbit
            self.draw_orbit(orbit_radius)

            # Position planet (for now, just place them at different positions)
            angle = math.radians(orbit_radius * 30)  # Different angles for each planet
            x = orbit_radius * math.cos(angle)
            z = orbit_radius * math.sin(angle)

            glTranslatef(x, 0, z)
            glColor3f(color[0], color[1], color[2])
            gluSphere(gluNewQuadric(), planet_size, 15, 15)
            self.draw_planet_center(0, 0, 0)

            glPopMatrix()

    def draw_text_overlay(self) -> None:
        """Draw text overlay using pygame surface rendered to OpenGL texture"""
        # Create a pygame surface for text rendering
        text_surface = pygame.Surface((350, 250), pygame.SRCALPHA)
        text_surface.fill((0, 0, 0, 180))  # Semi-transparent black background

        # Render text to the surface
        instructions = [
            "3D Camera Controls:",
            "",
            "SPACE - Toggle mouse capture",
            "WASD - Move horizontally",
            "Q/E - Move up/down",
            "Mouse - Look around",
            "Mouse wheel - Zoom",
            "ESC - Exit",
            "",
            f"Mouse: {'CAPTURED' if self.camera.mouse_captured else 'FREE'}"
        ]

        y_offset = 10
        line_height = 22

        for i, line in enumerate(instructions):
            if line.strip():  # Don't render empty lines
                if line.startswith("3D Camera"):
                    color = (255, 255, 100)  # Yellow for title
                elif line.startswith("Mouse:"):
                    color = (100, 255, 100) if self.camera.mouse_captured else (255, 100, 100)
                else:
                    color = (255, 255, 255)  # White for instructions

                text_img = self.font.render(line, True, color)
                text_surface.blit(text_img, (10, y_offset + i * line_height))

        # Convert pygame surface to OpenGL texture - flip vertically for OpenGL
        text_data = pygame.image.tostring(pygame.transform.flip(text_surface, False, True), "RGBA", False)

        # Save OpenGL state
        glMatrixMode(GL_PROJECTION)
        glPushMatrix()
        glLoadIdentity()
        glOrtho(0, self.width, 0, self.height, -1, 1)

        glMatrixMode(GL_MODELVIEW)
        glPushMatrix()
        glLoadIdentity()

        # Disable depth test and lighting for 2D overlay
        glDisable(GL_DEPTH_TEST)
        glDisable(GL_LIGHTING)

        # Enable blending and texturing
        glEnable(GL_BLEND)
        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA)
        glEnable(GL_TEXTURE_2D)

        # Create and bind texture
        texture_id = glGenTextures(1)
        glBindTexture(GL_TEXTURE_2D, texture_id)
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR)
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR)
        glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, text_surface.get_width(),
                     text_surface.get_height(), 0, GL_RGBA, GL_UNSIGNED_BYTE, text_data)

        # Draw textured quad - positioned correctly for top-left corner
        glColor4f(1.0, 1.0, 1.0, 1.0)
        glBegin(GL_QUADS)
        glTexCoord2f(0, 0)
        glVertex2f(10, self.height - 10 - text_surface.get_height())
        glTexCoord2f(1, 0)
        glVertex2f(10 + text_surface.get_width(), self.height - 10 - text_surface.get_height())
        glTexCoord2f(1, 1)
        glVertex2f(10 + text_surface.get_width(), self.height - 10)
        glTexCoord2f(0, 1)
        glVertex2f(10, self.height - 10)
        glEnd()

        # Clean up texture
        glDeleteTextures([texture_id])
        glDisable(GL_TEXTURE_2D)
        glDisable(GL_BLEND)

        # Restore OpenGL state
        glEnable(GL_DEPTH_TEST)
        glEnable(GL_LIGHTING)

        glPopMatrix()
        glMatrixMode(GL_PROJECTION)
        glPopMatrix()
        glMatrixMode(GL_MODELVIEW)

    def draw_coordinate_axes(self) -> None:
        """Draw coordinate axes for reference"""
        glDisable(GL_LIGHTING)
        glLineWidth(3)
        glBegin(GL_LINES)

        # X axis - Red
        glColor3f(1, 0, 0)
        glVertex3f(0, 0, 0)
        glVertex3f(3, 0, 0)

        # Y axis - Green
        glColor3f(0, 1, 0)
        glVertex3f(0, 0, 0)
        glVertex3f(0, 3, 0)

        # Z axis - Blue
        glColor3f(0, 0, 1)
        glVertex3f(0, 0, 0)
        glVertex3f(0, 0, 3)

        glEnd()
        glLineWidth(1)
        glEnable(GL_LIGHTING)

    def handle_events(self) -> bool:
        """Handle all pygame events"""
        for event in pygame.event.get():
            if event.type == pygame.QUIT:
                return False

            elif event.type == pygame.KEYDOWN:
                if event.key == pygame.K_ESCAPE:
                    return False
                elif event.key == pygame.K_SPACE:
                    self.camera.toggle_mouse_capture()

            elif event.type == pygame.MOUSEMOTION:
                self.camera.handle_mouse_motion(event.pos, pygame.mouse.get_pressed())

            elif event.type == pygame.MOUSEWHEEL:
                self.camera.handle_mouse_wheel(event.y)

        return True

    def update(self, dt: float) -> None:
        """Update game state"""
        if not isinstance(dt, (int, float)) or dt < 0:
            raise ValueError("dt must be a non-negative numeric value")
        keys = pygame.key.get_pressed()
        self.camera.handle_keyboard(keys, dt)

    def render(self) -> None:
        """Render the scene"""
        # Clear the screen
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)

        # Reset modelview matrix
        glLoadIdentity()

        # Apply camera transformation
        self.camera.apply_camera_transform()

        # Draw coordinate axes at origin
        # self.draw_coordinate_axes()

        # Render solar system in front-center, lower
        glPushMatrix()
        glTranslatef(0, -3, -6)  # Center, down, and forward from camera
        self.generate_solar_system()
        glPopMatrix()

        # Draw text overlay on top
        self.draw_text_overlay()

        pygame.display.flip()

    def run(self) -> None:
        """Main game loop"""
        running = True
        print("3D Camera Controls Loaded!")
        print("Controls:")
        print("  SPACE - Toggle mouse capture")
        print("  WASD - Move horizontally")
        print("  Q/E - Move up/down")
        print("  Mouse - Look around (when captured)")
        print("  Mouse wheel - Zoom")
        print("  ESC - Exit")
        print(f"\nCamera starts at {self.camera.position} looking toward origin")
        print("You should see cubes and a solar system!")

        while running:
            dt = self.clock.tick(60) / 1000.0  # Delta time in seconds

            running = self.handle_events()
            self.update(dt)
            self.render()

        pygame.quit()


if __name__ == "__main__":
    # Test Camera3D initialization first
    try:
        print("Testing Camera3D initialization...")
        test_camera = Camera3D()
        print("✓ Default Camera3D initialization successful")
        
        test_camera2 = Camera3D([1, 2, 3])
        print("✓ Custom position Camera3D initialization successful")
        
        test_camera3 = Camera3D([1, 2, 3], [0, 0, -1], [0, 1, 0])
        print("✓ Full custom Camera3D initialization successful")
        
    except Exception as e:
        print(f"✗ Camera3D Error: {e}")
        exit(1)
    
    # Test Renderer initialization
    try:
        print("Testing Renderer initialization...")
        renderer = Renderer()
        print("✓ Renderer initialization successful")
        renderer.run()
    except Exception as e:
        print(f"✗ Renderer Error: {e}")
        exit(1)